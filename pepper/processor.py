import os
import numpy as np
import awkward as ak
import uproot
import coffea.processor
from coffea.nanoevents import NanoAODSchema
import h5py
import json
import logging
from time import time
import abc
import uuid
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import concurrent.futures
from tqdm import tqdm

import pepper
from pepper import Selector, OutputFiller, HDF5File
import pepper.config
import pepper.htcondor


logger = logging.getLogger(__name__)


class Processor(coffea.processor.ProcessorABC):
    """Class implementing input/output, setup of histograms, and utility
    classes

    It implements many technicalities but no functionality related to physics.
    Classes deriving from it are supposed to implement cuts or particle
    definitions.

    Attributes
    ----------
    config_class
        Class to use for config parsing
    schema_class
        Class to use as schema for the input data (usually NanoAODSchema)
    """
    config_class = pepper.Config
    schema_class = NanoAODSchema

    def __init__(self, config, eventdir):
        """
        Parameters
        ----------
        config
            Instance of ``config_class``, containing the configuration to use
        eventdir
            Path to the destination directory, where the per event output is
            saved. Every chunk will be saved in its own file. If `None`,
            nothing will be saved.
        """
        self._check_config_integrity(config)
        self.config = config
        if eventdir is not None:
            self.eventdir = os.path.realpath(eventdir)
        else:
            self.eventdir = None

        self.rng_seed = self._load_rng_seed()
        self.loglevel = logging.getLogger("pepper").level

    @staticmethod
    def _check_config_integrity(config):
        """Is called when initialized and is supposed to check the
        configuration for obvious errors, so that the user has an immediate
        error message.
        """
        # Check for duplicate names in columns_to_save
        column_names = []
        if "columns_to_save" in config:
            to_save = config["columns_to_save"]
            if isinstance(to_save, dict):
                spec = to_save.items()
            else:
                spec = zip([None] * len(to_save), to_save)
            for key, specifier in spec:
                if key is None:
                    datapicker = pepper.hist_defns.DataPicker(specifier)
                    key = datapicker.name
                if key in column_names:
                    raise pepper.config.ConfigError(
                        f"Ambiguous column to save '{key}' (from {specifier})")
                else:
                    column_names.append(key)

    @staticmethod
    def _get_hists_from_config(config, key, todokey):
        """Get all histograms to create from config. The config allows the
        specification of a list of histograms to do, even if there are more
        histograms defined in the config."""
        if key in config:
            hists = config[key]
        else:
            hists = {}
        if todokey in config and len(config[todokey]) > 0:
            new_hists = {}
            for name in config[todokey]:
                if name in hists:
                    new_hists[name] = hists[name]
            hists = new_hists
            logger.info("Doing only the histograms: " +
                        ", ".join(hists.keys()))
        return hists

    def _load_rng_seed(self):
        """Load the random number generator seed. The seed is a large integer
        saved in a txt file. The location of the txt file is obtained from
        the configuration. If it does not exists, a new seed is made and saved
        to the txt file."""
        if "rng_seed_file" not in self.config:
            return np.random.SeedSequence().entropy
        seed_file = self.config["rng_seed_file"]
        if os.path.exists(seed_file):
            with open(seed_file) as f:
                try:
                    seed = int(f.read())
                except ValueError as e:
                    raise pepper.config.ConfigError(
                        f"Not an int in rng_seed_file '{seed_file}'")\
                        from e
                return seed
        else:
            rng_seed = np.random.SeedSequence().entropy
            with open(seed_file, "w") as f:
                f.write(str(rng_seed))
            return rng_seed

    def preprocess(self, datasets):
        """Modify the list of data sets that are processed

        The main purpose is method is when subclasses should be run only
        on specific data sets. These subclasses can enfored this here

        Parameters
        ----------
        datasets
            Dict mapping data set names to list of data set sources. Sources
            can be path or full CMS data set names.

        Returns
        -------
        datasets
            The same as the ```datasets``` parameters, but modified if needed
            by the processor.
        """
        return datasets

    @staticmethod
    def postprocess(accumulator):
        """Modify the output of the produced processor

        This could be overwritten by subclasses if they want to modify the
        output slightly

        Parameters
        ----------
        accumulator
            The output of the processor

        Returns
        -------
        accumulator
            Modified, if nessecary, version of the input ``accumulator``
        """
        return accumulator

    def _open_output(self, dsname, filetype):
        """Try to open an output file for writing the per-event data. It is
        ensured the file is newly created and does not overwrite existing data.

        Parameters
        ----------
        dsname
            Name of the data set of the data
        filetype
            Either "root" or "hdf5". The type of the output file

        Returns
        -------
            File object of the opened file
        """
        if filetype == "root":
            ext = ".root"
        elif filetype == "hdf5":
            ext = ".h5"
        else:
            raise ValueError(f"Invalid filetype: {filetype}")
        dsname = dsname.replace("/", "_")
        dsdir = os.path.join(self.eventdir, dsname)
        os.makedirs(dsdir, exist_ok=True)
        i = 0
        while True:
            filepath = os.path.join(dsdir, str(i).zfill(4) + ext)
            try:
                f = open(filepath, "x")
            except (FileExistsError, OSError):
                i += 1
                continue
            else:
                break
        f.close()
        if filetype == "root":
            f = uproot.recreate(filepath)
        elif filetype == "hdf5":
            f = h5py.File(filepath, "w")
        logger.debug(f"Opened output {filepath}")
        return f

    def _prepare_saved_columns(self, selector):
        """Creates an array to be saved as per-event data. The content is taken
        from the selectors and the data pickers defined in the config."""
        columns = {}
        if "columns_to_save" in self.config:
            to_save = self.config["columns_to_save"]
        else:
            to_save = []
        if isinstance(to_save, dict):
            spec = to_save.items()
        else:
            spec = zip([None] * len(to_save), to_save)
        for key, specifier in spec:
            datapicker = pepper.hist_defns.DataPicker(specifier)
            item = datapicker(selector.data)
            if item is None:
                logger.info("Skipping to save column because it is not "
                            f"present: {specifier}")
                continue
            if key is None:
                key = datapicker.name
            columns[key] = item
        return ak.Array(columns)

    def _prepare_saved_categories(self, selector):
        cat_dict = defaultdict(dict)
        for cat_name, cats in selector.cats.items():
            for cat in cats:
                cat_dict[cat_name][cat] = selector.data[cat]
        return ak.Array(cat_dict)

    def _save_per_event_info_hdf5(
            self, dsname, selector, identifier, save_full_sys=True):
        """Save the per-event info into an HDF5 file"""
        out_dict = {"dsname": dsname, "identifier": identifier}
        out_dict["events"] = self._prepare_saved_columns(selector)
        cutnames, cutflags = selector.get_cuts()
        out_dict["cutnames"] = cutnames
        out_dict["cutflags"] = cutflags
        if ("save_categories_per_event" not in self.config or
                self.config["save_categories_per_event"]):
            out_dict["categories"] = \
                self._prepare_saved_categories(selector)
        if (self.config["compute_systematics"] and save_full_sys
                and selector.systematics is not None):
            out_dict["systematics"] = ak.flatten(selector.systematics,
                                                 axis=0)
        elif selector.systematics is not None:
            out_dict["weight"] = ak.flatten(
                selector.systematics["weight"], axis=0)
        with self._open_output(dsname, "hdf5") as f:
            outf = HDF5File(f)
            for key in out_dict.keys():
                outf[key] = out_dict[key]

    @staticmethod
    def _separate_masks_for_root(arrays):
        """Seperate a masked awkward array into an unmasked array and an array
        defining its mask.

        Parameters
        ----------
        arrays
            Dict of awkard arrays to unmask

        Returns
        -------
            Dict with unsmaked arrays and their masks. The masks have the same
            key prefixed by "mask"
        """
        ret = {}
        for key, array in arrays.items():
            if (not isinstance(array, ak.Array)
                    or not pepper.misc.akismasked(array)):
                ret[key] = array
                continue
            if array.ndim > 2:
                raise ValueError(
                    f"Array '{key}' as too many dimensions for ROOT output")
            if "mask" + key in arrays:
                raise RuntimeError(f"Output named 'mask{key}' already present "
                                   "but need this key for storing the mask")
            ret["mask" + key] = ~ak.is_none(array)
            if len(array.fields) > 0:
                ret[key] = ak.fill_none(array, {k: 0 for k in array.fields})
            else:
                ret[key] = ak.fill_none(array, 0)
        return ret

    def _save_per_event_info_root(self, dsname, selector, identifier,
                                  save_full_sys=True):
        """Save the per-event info into a Root file"""
        out_dict = {"dsname": dsname, "identifier": str(identifier)}
        cutnames, cutflags = selector.get_cuts()
        out_dict["Cutnames"] = str(cutnames)
        if len(selector.data) > 0:
            events = self._prepare_saved_columns(selector)
            # Workaround: Use ak.packed to make sure offset arrays of virtual
            # arrays are not given to uproot. Uproot has a bug for these.
            events = {f: ak.packed(events[f]) for f in ak.fields(events)}
            additional = {}
            additional["cutflags"] = cutflags
            if selector.systematics is not None:
                additional["weight"] = selector.systematics["weight"]
                if self.config["compute_systematics"] and save_full_sys:
                    for field in ak.fields(selector.systematics):
                        additional[f"systematics_{field}"] = \
                            selector.systematics[field]

            for key in additional.keys():
                if key in events:
                    raise RuntimeError(
                        f"branch named '{key}' already present in Events tree")
            events.update(additional)
            events = self._separate_masks_for_root(events)
            out_dict["Events"] = events
            if ("save_categories_per_event" not in self.config or
                    self.config["save_categories_per_event"]):
                cats = self._prepare_saved_categories(selector)
                for cat in ak.fields(cats):
                    out_dict[f"Categories/{cat}"] = \
                        self._separate_masks_for_root(
                            {f: ak.packed(cats[cat][f])
                             for f in ak.fields(cats[cat])})
        with self._open_output(dsname, "root") as outf:
            for key in out_dict.keys():
                outf[key] = out_dict[key]

    def save_per_event_info(self, dsname, selector, save_full_sys=True):
        """Save the per-event info

        Parameters
        ----------
        dsname
            Name of the data set of the data
        selector
            Selector containing data, systematics and all the other info
            we save
        identifier
            Touple that uniquely identifies the data that goes into the file
        save_full_sys
            Whether to save all systematic variations. If ``False`` only
            the event weight is saved
        """
        idn = self.get_identifier(selector)
        logger.debug("Saving per event info")
        if "column_output_format" in self.config:
            outformat = self.config["column_output_format"].lower()
        else:
            outformat = "root"
        if outformat == "root":
            self._save_per_event_info_root(
                dsname, selector, idn, save_full_sys)
        elif outformat == "hdf5":
            self._save_per_event_info_hdf5(
                dsname, selector, idn, save_full_sys)
        else:
            raise pepper.config.ConfigError(
                "Invalid value for column_output_format, must be 'root' "
                "or 'hdf'")

    @staticmethod
    def get_identifier(data):
        """Get a unique identifier for the data as used in the per-event data
        file

        Parameters
        ----------
        data
            Data array (usually NanoEvents) or Selector

        Returns
        -------
            Tuple uniquely identifing the data
        """
        meta = data.metadata
        return meta["filename"], meta["entrystart"], meta["entrystop"]

    def process(self, data):
        """Do all setup steps of the selector, output filler, follwed by
        performing the actual selection and saving the output

        Parameters
        ----------
        data
            Data array (usually NanoEvents)

        Returns
        -------
            Output from the processor, containing hists and/or cutflows
        """
        pepper_logger = logging.getLogger("pepper")
        try:
            jobad = pepper.htcondor.get_htcondor_jobad()
        except OSError:
            pass
        else:
            if not getattr(pepper_logger, "is_on_condor", False):
                pepper_logger.addHandler(logging.StreamHandler())
                pepper_logger.setLevel(self.loglevel)
                pepper_logger.is_on_condor = True

                jname, jid, jtime = jobad["GlobalJobId"].split("#", 2)
                logger.debug(f"Running on machine {jobad['RemoteHost']} for "
                             f"job {jid}")

        try:
            return self._process_inner(data)
        except Exception as e:
            if getattr(pepper_logger, "is_on_condor", False):
                logger.exception(e)
            raise

    def _process_inner(self, data):
        """Inner part of the ``process()`` method, so that it can easily be
        part of a try-block"""
        starttime = time()
        dsname = data.metadata["dataset"]
        filename = data.metadata["filename"]
        entrystart = data.metadata["entrystart"]
        entrystop = data.metadata["entrystop"]
        logger.debug(f"Started processing {filename} from event "
                     f"{entrystart} to {entrystop - 1} for dataset {dsname}")
        is_mc = (dsname in self.config["mc_datasets"].keys())

        filler = self.setup_outputfiller(dsname, is_mc)
        selector = self.setup_selection(data, dsname, is_mc, filler)
        self.process_selection(selector, dsname, is_mc, filler)

        if self.eventdir is not None:
            self.save_per_event_info(dsname, selector)

        timetaken = time() - starttime
        logger.debug(f"Processing finished. Took {timetaken:.3f} s.")
        return filler.output

    def setup_outputfiller(self, dsname, is_mc):
        """Create a new output filler to be used throughout the selection. The
        output filler is responsible to create the output of the processor,
        including histograms and cutflows

        Parameters
        ----------
        dsname
            Name of the data set that is processed
        is_mc
            Whether the data is simulation

        Returns
        -------
            An instance of an ``OutputFiller`` to be used for the selection
        """
        sys_enabled = self.config["compute_systematics"]

        if dsname in self.config["dataset_for_systematics"]:
            dsname_in_hist = self.config["dataset_for_systematics"][dsname][0]
            sys_overwrite = self.config["dataset_for_systematics"][dsname][1]
        elif ("datasets_to_group" in self.config
              and dsname in self.config["datasets_to_group"]):
            dsname_in_hist = self.config["datasets_to_group"][dsname]
            sys_overwrite = None
        else:
            dsname_in_hist = dsname
            sys_overwrite = None

        if "cuts_to_histogram" in self.config:
            cuts_to_histogram = self.config["cuts_to_histogram"]
        else:
            cuts_to_histogram = None

        hists = self._get_hists_from_config(
            self.config, "hists", "hists_to_do")
        filler = OutputFiller(
            hists, is_mc, dsname, dsname_in_hist, sys_enabled,
            sys_overwrite=sys_overwrite, cuts_to_histogram=cuts_to_histogram)

        return filler

    def setup_selection(self, data, dsname, is_mc, filler):
        """Create a new selector that is to be used throughout the selection.
        The selector lets us specify cuts and new columns.

        Parameters
        ----------
        data
            Data array (usually NanoEvents)
        dsname
            Name of the data set that is processed
        is_mc
            Whether the data is simulation
        filler
            The output filler used in the selection

        Returns
        -------
            A new instance of ``Selector`` for the selection
        """
        if is_mc:
            if (("norm_genweights" in self.config
                    and self.config["norm_genweights"])
                    or ("genweights_to_norm" in self.config and
                        dsname in self.config["genweights_to_norm"])):
                genweight = np.sign(data["genWeight"])
            else:
                genweight = data["genWeight"]
        else:
            genweight = None
        # Use a different seed for every chunk in a reproducable way
        seed = (self.rng_seed, uuid.UUID(data.metadata["fileuuid"]).int,
                data.metadata["entrystart"])
        selector = Selector(data, genweight, filler.get_callbacks(),
                            rng_seed=seed)
        return selector

    @abc.abstractmethod
    def process_selection(self, selector, dsname, is_mc, filler):
        """Do selection steps, e.g. cutting, defining objects.

        This is to be defined in the practial implementations of the
        processors. Users that want to implement an analysis should inherit in
        some way from the processor and overwrite this method.

        Parameters
        ----------
        selector
            A pepper.Selector object with the event data
        dsname
            Name of the data set that is processed
        is_mc
            Whether the data is simulation
        filler
            pepper.OutputFiller object to controll how the output is structured
        """

    @staticmethod
    def _get_cuts(output):
        """Get a list of cuts in the order they are applied

        The cuts are obtained from an output's cutflow.

        Paramters
        ---------
        output
            The output in which the cutflow is found

        Returns
        -------
            List of cuts

        Raises
        ------
        ValueError
            When no ordering of cuts could be identified
        """
        cutflow_all = output["cutflows"]
        cut_lists = [list(cutflow.keys()) for cutflow
                     in cutflow_all.values()]
        cuts_precursors = defaultdict(set)
        for cut_list in cut_lists:
            for i, cut in enumerate(cut_list):
                cuts_precursors[cut].update(set(cut_list[:i]))
        cuts = []
        while len(cuts_precursors) > 0:
            for cut, precursors in cuts_precursors.items():
                if len(precursors) == 0:
                    cuts.append(cut)
                    for p in cuts_precursors.values():
                        p.discard(cut)
                    cuts_precursors.pop(cut)
                    break
            else:
                raise ValueError("No well-defined ordering of cuts "
                                 "for all datasets found")
        return cuts

    @staticmethod
    def _prepare_cutflows(proc_output):
        """Convert the cutflows into a dictionary. Cutflows are produced as
        one bin histograms, thus conversion is needed. Aditionally, this adds
        a sum (with the key "all")."""
        cutflows = proc_output["cutflows"]
        output = {}
        for dataset, cf1 in cutflows.items():
            output[dataset] = {"all": defaultdict(float)}
            for cut, cf2 in cf1.items():
                cf = pepper.misc.get_hist_cat_values(cf2)
                for cat_position, value in cf.items():
                    output_for_cat = output[dataset]
                    for cat_coordinate in cat_position:
                        if cat_coordinate not in output_for_cat:
                            output_for_cat[cat_coordinate] = {}
                        output_for_cat = output_for_cat[cat_coordinate]
                    if len(cat_position) > 0:
                        output_for_cat[cut] = value.sum()
                    output[dataset]["all"][cut] += value.sum()
        return output

    @staticmethod
    def _save_histograms_inner(key, histdict, cuts, hist_col, format):
        """Save a histogram

        This method does the actual work and can be run
        in paralel. This may take some time due to having to sum histograms
        across different data sets and in case of the Root format, having
        to split into sub-histograms.

        Parameters
        ----------
        key
            Key to be used in ``hist_col``. Tuple with the first element being
            the cut name
        histdict
            The histogram split into sub-histograms, one for each data set
        cuts
            List of cuts that have been applied
        hist_col
            HistCollection instance, which is used to save the histogram
        format
            Either "hist" or "root". The format to save the histogram in.
            Usually "hist" is faster.

        Returns
        -------
            The ``hist_col``

        """
        hist_sum = None
        cats_present = set()
        for dataset, hist in histdict.items():
            cats_present |= hist_col.get_cats_present(hist)
            if hist_sum is None:
                hist_sum = hist.copy()
            else:
                hist_sum += hist
        cutnum = cuts.index(key[0])
        if format == "root":
            ext = ".root"
        elif format == "hist":
            ext = ".coffea"
        else:
            raise ValueError(f"Invalid hist format: {format}")
        fname = f"Cut_{cutnum:03}_{'_'.join(key)}{ext}"
        fname = fname.replace("/", "")
        key = key + (None,) * (len(hist_col.key_fields) - len(key))
        hist_col.save(key, hist_sum, fname, format, cats_present=cats_present)
        return hist_col

    @classmethod
    def save_histograms(cls, format, output, dest, threads=10):
        """Save histograms to files

        Parameters
        ----------
        format
            Either "hist" or "root". The format to save the histogram in.
            Usually "hist" is faster.
        output
            Output from the processor's output filler
        dest
            Path to the destination directory to save the histograms in
        threads
            Number of processes to run in parallel to do the saving
        """
        cuts = cls._get_cuts(output)
        hists = defaultdict(dict)
        data = {"cuts": cuts}
        hist_col = pepper.HistCollection(dest, ["cut", "hist"], userdata=data)
        for dataset, hists_per_ds in output["hists"].items():
            for key, hist in hists_per_ds.items():
                hists[key][dataset] = hist
        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = []
            for key, histdict in hists.items():
                futures.append(executor.submit(
                    cls._save_histograms_inner, key, histdict, cuts, hist_col,
                    format))

            for future in tqdm(concurrent.futures.as_completed(futures),
                               desc="Saving histograms", total=len(futures)):
                hist_col += future.result()
        with open(os.path.join(dest, "hists.json"), "w") as f:
            hist_col.save_metadata_json(f)

    def save_output(self, output, dest):
        """Save the histograms and cutflows to files

        Parameters
        ----------
        output
            Output from the processor's output filler
        dest
            Destination direction
        """
        # Save cutflows
        with open(os.path.join(dest, "cutflows.json"), "w") as f:
            json.dump(self._prepare_cutflows(output), f, indent=4)

        if "histogram_format" in self.config:
            hform = self.config["histogram_format"].lower()
        else:
            hform = "hist"
        if hform not in ["coffea", "root", "hist"]:
            logger.warning(
                f"Invalid histogram format: {hform}. Saving as hist")
            hform = "hist"
        hist_dest = os.path.join(dest, "hists")
        os.makedirs(hist_dest, exist_ok=True)
        self.save_histograms(hform, output, hist_dest)
