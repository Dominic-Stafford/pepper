import os
import abc
from dataclasses import dataclass, fields
from typing import Optional
import threading
import queue
from copy import deepcopy
from functools import partial
import time
import logging
import uproot
import coffea.processor
import coffea.util
from coffea.processor import set_accumulator
from coffea.processor.executor import (
    _compression_wrapper, _decompress, UprootMissTreeError, FileMeta)
from coffea.processor.accumulator import iadd as accum_iadd
from coffea.processor import ProcessorABC
from coffea.nanoevents import NanoEventsFactory
import cloudpickle
import uuid
import lz4.frame as lz4f
from tqdm import tqdm
import pepper


STATEFILE_VERSION = 3

logger = logging.getLogger(__name__)


class StateFileError(Exception):
    pass


def _wrap_execution(function, item):
    return item, function(item)


@dataclass
class ResumableExecutor(abc.ABC, coffea.processor.executor.ExecutorBase):
    """Abstract base class for executors that save their state and thus are
    able to resume if there was an interruption

    Attributes
    ----------
    state_file_name
        Name of the file to save and read the state to/from
    remove_state_at_end
        Bool, if true, remove the state file after successful completion
    save_interval
        Seconds that have to pass before the state is saved after start or the
        last save. The sate is saved only after the completion of an item
    """

    state_file_name: Optional[str] = None
    remove_state_at_end: bool = False
    save_interval: int = 60  # seconds

    def __post_init__(self):
        self.state = {"items_done": [], "accumulator": None,
                      "version": STATEFILE_VERSION, "userdata": {}}
        # Thread to manage accumulation and state saving
        self._state_manager = threading.Thread(target=self._manage_state)
        # To communicate with the thread
        self._state_manager_queue = queue.Queue()
        self._is_running = threading.Event()
        self._has_exception = threading.Event()
        # tqdm progress bar
        self._progress = None

    def _manage_state(self):
        """Method run inside a separate thread, which watches the
        ``_state_manager_queue`` for new results, adds them to the ones already
        in the state and saves the state to a file if necessary
        """
        state = self.state
        state_changed = False
        nextstatebackup = time.time() + self.save_interval
        while not self._has_exception.is_set():
            try:
                result = self._state_manager_queue.get(timeout=0.1)
            except queue.Empty:
                if not self._is_running.is_set():
                    break
            else:
                state["items_done"], state["accumulator"] = self._accumulate(
                    [result], state["items_done"], state["accumulator"])
                state_changed = True
                self._state_manager_queue.task_done()
                if self._progress is not None:
                    self._progress.update(1)
            if nextstatebackup <= time.time() and state_changed:
                self.save_state()
                state_changed = False
                nextstatebackup = time.time() + self.save_interval
        if state_changed:
            self.save_state()

    def _accumulate(self, results, items_done=None, accumulator=None):
        """Accumulates results and adds them to an accumulator

        Accumulation essentially uses the += operator on the objects inside the
        results.

        Parameters
        ----------
        results
            Results as they come from the workers. Will be decompressed if
            neccessary.
        items_done
            The items in results will be added to this list
        accumulator
            The results will be added to this accumulator

        Returns
        -------
        items_done
            Same as the ``items_done`` parameter or a new list if it was
            ``None``
        accumulator
            Same as the ``accumulator`` parameter or a new accumulator if it
            was ``None``
        """
        if items_done is None:
            items_done = []
        for result in results:
            if self.compression is not None:
                result = _decompress(result)
            item, res = result
            items_done.append(item)
            if accumulator is None:
                accumulator = res
            else:
                accum_iadd(accumulator, res)
        return items_done, accumulator

    def copy(self, **kwargs):
        """Create a shallow copy of the executor

        This is similar to Coffea's ExecutorBase.copy.

        Parameters
        ----------
        **kwargs
            Attributes that will be overwritten in the copy

        Returns
        -------
            A shallow copy of `self`
        """
        tmp = {f.name: getattr(self, f.name) for f in fields(self)}
        tmp.update(kwargs)
        instance = type(self)(**tmp)
        # Need deep copy here to not modify the accumulator later
        instance.state = deepcopy(self.state)
        return instance

    def load_state(self, filename=None):
        """Load a previous state from a file

        Parameters
        ----------
        filename
            Name of the file to read the state from. If not given, will load
            from ``self.state_file_name``.
        """

        if filename is None:
            if self.state_file_name is None:
                raise ValueError("Need filename")
            filename = self.state_file_name
        state = coffea.util.load(filename)
        if state["version"] != STATEFILE_VERSION:
            raise StateFileError("State file made by incompatible version: "
                                 f"{filename}")
        self.state = state

    def reset_state(self):
        """Set the state to an empty one
        """
        self.state = {"items_done": [], "accumulator": None,
                      "version": STATEFILE_VERSION, "userdata": {}}

    def __call__(self, items, function, accumulator):
        """Start the execution, usually to be called from inside the Coffea
        Runner

        Parameters
        ----------
        items
            Each item is hashable and uniquely identifies each task to be done
        function
            The function to be executied on ``items``
        accumulator
            Accumulator to use. Can contain results already or can also be
            ``None``
        """
        items_done = set(self.state["items_done"])
        items = [item for item in items if item not in items_done]
        if len(items) == 0:
            return self.state["accumulator"], 0

        if accumulator is not None and self.state["accumulator"] is not None:
            accumulator.add(self.state["accumulator"])
            self.state["accumulator"] = accumulator
        elif accumulator is not None:
            self.state["accumulator"] = accumulator
        elif self.state["accumulator"] is not None:
            accumulator = self.state["accumulator"]

        function = partial(_wrap_execution, function)
        if self.compression is not None:
            function = _compression_wrapper(self.compression, function)

        gen = self._submit(items, function)
        with tqdm(
            total=len(items),
            desc=self.desc,
            unit=self.unit,
            disable=not self.status
        ) as self._progress:
            self._is_running.set()
            self._state_manager.start()
            try:
                while True:
                    result = next(gen)
                    self._state_manager_queue.put(result)
            except StopIteration:
                pass
            except (Exception, KeyboardInterrupt):
                self._has_exception.set()
                raise
            finally:
                gen.close()
            self._is_running.clear()
            self._state_manager.join()
        if (self.state_file_name is not None
                and self.remove_state_at_end
                and os.path.exists(self.state_file_name)):
            os.remove(self.state_file_name)
        res = self.state["accumulator"]
        return res, 0

    @abc.abstractmethod
    def _submit(self, items, function):
        """Defines how to submit and execute the work

        Parameters
        ----------
        items
            Each item is so be called on ``function``
        function
            Function to be executed

        Returns
        -------
        results
            Generator for the results. Results that are completed first should
            also be yielded first
        """
        return

    def save_state(self):
        """Save the state of this executor to the disk

        The file created is named according to the ``state_file_name``
        attribute and cam be used to obtain an executor of the same state using
        ``load_state``.
        """
        # Save state to a new file and only replace the previous state file
        # when writing is finished. This avoids leaving only an invalid state
        # file if the program is terminated during writing.
        i = 0
        while True:
            output = os.path.join(os.path.dirname(self.state_file_name),
                                  f"pepper_temp{i}.coffea")
            if not os.path.exists(output):
                break
            i += 1
        coffea.util.save(self.state, output)
        os.replace(output, self.state_file_name)


# This data class allows us to put the cluster parameter first in
# ClusterExecutor.__init__
@dataclass
class _WithCluster:
    """
    Attributes
    ----------
    cluster
        To be used to execute work
    """
    cluster: pepper.htcondor.Cluster


@dataclass
class ClusterExecutor(ResumableExecutor, _WithCluster):
    """Uses ``pepper.htcondor.Cluster`` to execute work
    """
    @staticmethod
    def get_taskname(item, i):
        if hasattr(item, "entrystart"):
            return (f"{item.dataset[:40]}/{os.path.basename(item.filename)}/"
                    f"{item.entrystart}:{item.entrystop}/chunk{i:06}")
        else:
            return (f"{item.dataset[:40]}/{os.path.basename(item.filename)}/"
                    f"/chunk{i:06}")

    def _submit(self, items, function):
        tasknames = list(map(self.get_taskname, items, range(len(items))))
        yield from self.cluster.process(function, items, key=tasknames)


class Runner(coffea.processor.Runner):
    """
    This is used to make it possible retry a different XRootD server if the
    first server gave an error

    Making a subclass is a bad solution but there seems to be no other solution
    other than making our own Runner. Coffea's Runner is using ``uproot.open``
    inside ``metadata_fetcher`` and ``_work_function``.
    """
    @staticmethod
    def resolve_lfn(lfn, store, xrootddomain, skippaths, url_blacklist=None):
        """Converts logical file names (LFNs) to physical file names that can
        be understood by ``uproot.open``

        Parameters
        ----------
        lfn
            If it starts with 'cmslfn://', it is interpreseted as logical file
            name, otherwise it is assumed it already is a physical file name
        store
            Path to the store directory for local access to the file
        xrootddomain
            Domain of the redirector server to find sites that offer the file
            via XRootD
        skippaths
            Blacklist of physical file paths to ignore
        url_blacklist
            Optional; a blacklist of XRootD URLs which should not be used for
            resolving the file

        Returns
        -------
        filepaths
            Physical file paths associated to ``lfn``
        """
        if lfn.startswith("cmslfn://"):
            filepaths = pepper.datasets.resolve_lfn(lfn, store, xrootddomain,
                                                    url_blacklist)
        else:
            filepaths = [lfn]
        if skippaths is not None:
            filepaths = [p for p in filepaths if p not in skippaths]

        return filepaths

    @staticmethod
    def metadata_fetcher(xrootdtimeout, align_clusters, item):
        filepaths = Runner.resolve_lfn(
            item.filename, item.metadata["store_path"],
            item.metadata["xrootddomain"], item.metadata["skippaths"],
            item.metadata["url_blacklist"])
        for filepath in filepaths:
            try:
                with uproot.open(
                        {filepath: None}, timeout=xrootdtimeout) as file:
                    try:
                        tree = file[item.treename]
                    except uproot.exceptions.KeyInFileError as e:
                        raise UprootMissTreeError(str(e)) from e

                    metadata = {}
                    if item.metadata:
                        metadata.update(item.metadata)
                    metadata.update({
                        "numentries": tree.num_entries,
                        "uuid": file.file.fUUID})
                    if align_clusters:
                        metadata["clusters"] = tree.common_entry_offsets()
                    out = set_accumulator(
                        [FileMeta(
                            item.dataset, item.filename, item.treename,
                            metadata)]
                    )
            except OSError as e:
                logger.warning(
                    "Got error while opening, continuing with alternative "
                    f"file location: {e}")
                continue
            break
        else:
            raise OSError(
                "None of the file paths found for the following file could be "
                f"opened: {item.filename}")

        return out

    @staticmethod
    def _work_function(
        format,
        xrootdtimeout,
        mmap,
        schema,
        cache_function,
        use_dataframes,
        savemetrics,
        item,
        processor_instance,
    ):
        if not isinstance(processor_instance, ProcessorABC):
            processor_instance = cloudpickle.loads(
                lz4f.decompress(processor_instance))
        # The ResumableExecutor might have loaded and old state, thus giving
        # old metadata possibly before the user changed the config.
        # Instead obtain the metadata from the processor_instance
        metadata = processor_instance.pepperitemmetadata

        filepaths = Runner.resolve_lfn(
            item.filename, metadata["store_path"],
            metadata["xrootddomain"], metadata["skippaths"],
            metadata["url_blacklist"])

        for filepath in filepaths:
            try:
                filecontext = uproot.open(
                    {filepath: None},
                    timeout=xrootdtimeout,
                    file_handler=uproot.MemmapSource
                    if mmap
                    else uproot.MultithreadedFileSource,
                )
            except OSError as e:
                logger.warning(
                    "Got error while opening, continuing with alternative "
                    f"file location: {e}")
                continue
            break
        else:
            raise OSError(
                "None of the file paths found for the following file could be "
                f"opened: {item.filename}")

        metadata = {
            "dataset": item.dataset,
            "filename": filepath,
            "treename": item.treename,
            "entrystart": item.entrystart,
            "entrystop": item.entrystop,
            "fileuuid": str(uuid.UUID(bytes=item.fileuuid))
            if len(item.fileuuid) > 0
            else "",
        }
        if item.usermeta is not None:
            metadata.update(item.usermeta)

        materialized = []
        with filecontext as file:
            factory = NanoEventsFactory.from_root(
                file=file,
                treepath=item.treename,
                entry_start=item.entrystart,
                entry_stop=item.entrystop,
                persistent_cache=cache_function(),
                schemaclass=schema,
                metadata=metadata,
                access_log=materialized,
            )
            events = factory.events()

            out = processor_instance.process(events)

        return {"out": out}
