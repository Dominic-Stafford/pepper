import os
import logging
import resource
import asyncio
import uuid
import dask
import dask.distributed
import dask_jobqueue
import shlex
from packaging.version import parse as parse_version


logger = logging.getLogger(__name__)


# Older versions of dask_jobqueue are missing the PR
# https://github.com/dask/dask-jobqueue/pull/610
# Workaround: Patch in necessary changes if version too old
if parse_version(dask_jobqueue.__version__) < parse_version("0.8.3"):
    from contextlib import suppress

    class HTCondorJob(dask_jobqueue.htcondor.HTCondorJob):
        async def _submit_job(self, script_filename):
            return await self._call(
                shlex.split(self.submit_command) + [script_filename])

        async def start(self):
            """Start workers and point them to our local scheduler"""

            with self.job_file() as fn:
                out = await self._submit_job(fn)
                self.job_id = self._job_id_from_submit_output(out)

            await super(dask_jobqueue.core.Job, self).start()

        async def close(self):
            await self._close_job(self.job_id, self.cancel_command)

        @classmethod
        async def _close_job(cls, job_id, cancel_command):
            if job_id:
                with suppress(RuntimeError):
                    await cls._call(shlex.split(cancel_command) + [job_id])

        @staticmethod
        async def _call(cmd, **kwargs):
            cmd_str = " ".join(cmd)

            proc = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
                **kwargs
            )

            out, err = await proc.communicate()
            out, err = out.decode(), err.decode()

            if proc.returncode != 0:
                raise RuntimeError(
                    "Command exited with non-zero exit code.\n"
                    "Exit code: {}\n"
                    "Command:\n{}\n"
                    "stdout:\n{}\n"
                    "stderr:\n{}\n".format(proc.returncode, cmd_str, out, err)
                )
            return out

    class HTCondorCluster(dask_jobqueue.htcondor.HTCondorCluster):
        job_cls = HTCondorJob
else:
    HTCondorCluster = dask_jobqueue.htcondor.HTCondorCluster


def get_site():
    """
    Returns
    -------
    hostname
        Name of the computing site currently on. If the site is
        unknown, returns the hostname
    """
    hostname = os.uname().nodename
    if hostname.endswith(".cern.ch") and hostname.startswith("lxplus"):
        return "lxplus"
    elif hostname.endswith(".desy.de") and hostname.startswith("naf-"):
        return "naf"
    else:
        return hostname


def get_dask_cluster(num_jobs, runtime=3*60*60, memory="2 GB", disk="3 GB",
                     cores=1, *, condorsubmit=None, condorenv=None,
                     logdir=None, memorylimit=0):
    """Get a Dask Jobqueue HTCondor cluster for a host

    Parameters
    ---------
    num_jobs
        Number of jobs to run in parallel
    runtime
        Requested runtime in seconds. If None, do not request a runtime
    memory
        Request memory. String with a unit like "GB" or an int.
    disk
        Request disk space. String with a unit like "GB" or an int.
    cores
        Total number of cores per job
    condorsubmit
        String containing additional parameters for the HTCondor submit file.
    condorenv
        Path to a Shell script that is sourced before the job starts
        If None, try to use the file pointed at by the local
        environment variable PEPPER_CONDOR_ENV and its contents
        instead. If PEPPER_CONDOR_ENV is also not set, no futher
        environment will be set up.
    logdir
        Directory where to store stdout and stderr logs
    memorylimit
        Maximum memory the job is allowed to use before it is killed by Dask.
        Default is 0, which sets no memory limit in Dask (leaving all memory
        management to Condor).

    Returns
    -------
    cluster
        Can be used within a distrubted Client to submit jobs to HTCondor
        Use ``client = distrubted.Client(cluster)``
    """

    site_config = {
        "lxplus": {
            # lxplus has a firewall in place, only allowing specific ports
            # https://batchdocs.web.cern.ch/specialpayload/dask.html
            "scheduler_options": {"port": 8786, "dashboard_address": ":0"},
            "worker_extra_args": ["--worker-port", "10000:10100"],
        }
    }

    # By default Dask Jobqueue sets the memory limit, at which Dask kills the
    # worker, to the same value as requested in the HTCondor submit file.
    # We don't want this, because usually one can use more memory than
    # requested through Condor. Thus explicitly set RequestMemory.
    memory = dask.utils.parse_bytes(memory)
    job_extra_directives = {
        "RequestMemory": str(int(memory / 2**20)),
        '+RequestRuntime': str(int(runtime))
    }
    if condorsubmit is not None:
        for param in condorsubmit.split("\n"):
            # Need to parse, HTCondorCluster only takes a dict
            key, val = param.split("=", 1)
            key = key.strip()
            val = val.strip()
            job_extra_directives[key] = val
    if condorenv is None and "PEPPER_CONDOR_ENV" in os.environ:
        condorenv = os.environ["PEPPER_CONDOR_ENV"]
    config = dict(
        name="PepperJob",
        cores=cores,
        memory=memorylimit,
        disk=disk,
        log_directory=logdir,
        job_extra_directives=job_extra_directives,
        job_script_prologue=["source " + shlex.quote(condorenv)],
        # Set port parameters to 0 to use random ports
        scheduler_options={"dashboard_address": ":0"}
    )
    site = get_site()
    config.update(site_config.get(site, {}))
    cluster = HTCondorCluster(**config)
    cluster.adapt(maximum_jobs=num_jobs)

    return cluster


def get_htcondor_jobad():
    """Get the HTCondor job AD as a dict of the job currently running in.
    If not running within a job, an OSError will be raised.
    For details on job AD see
    https://htcondor.readthedocs.io/en/latest/classad-attributes/job-classad-attributes.html
    """
    if "_CONDOR_JOB_AD" not in os.environ:
        raise OSError("Not inside HTCondor job")
    with open(os.environ["_CONDOR_JOB_AD"]) as f:
        jobad = f.readlines()
    ret = {}
    for line in jobad:
        k, v = line.split("=", 1)
        ret[k.strip()] = v.strip()
    return ret


class Cluster:
    """Run tasks either locally or on HTCondor.

    Internally, this class may call ``get_dask_cluster`` in order to submit
    jobs to HTCondor.
    """
    def __init__(
            self, num_jobs, condorsubmit=None, condorinit=None,
            logdir="pepper_logs", retries=None, condorsubmitfile=None,
            memory="2 GB", runtime=3*60*60):
        """
        Parameters
        ----------
        num_jobs
            The number of jobs to create on HTCondor. If ``None`` run locally
        condorsubmit
            Additional content to add to the HTCondor submit file
        condorinit
            Path to a script that will get sourced by the HTCondor jobs in
            order to initialize the environment
        logdir
            Directory to write log file to
        retries
            Number of retries if a job fails. If ``None`` retry indefinitely
        condorsubmitfile
            Path to a file containing additional content to add to the
            HTCondor submit file
        memory
            Request memory. String with a unit like "GB" or an int.
        """
        self.logdir = self.get_enumerated_dir(logdir)
        if num_jobs is None:
            # Run locally
            self.client = None
        else:
            # Run on HTCondor using Dask Jobqueue
            if condorsubmitfile is not None:
                if condorsubmit is None:
                    condorsubmit = ""
                with open(condorsubmitfile) as f:
                    condorsubmit += "\n" + f.read()
            dask_cluster = get_dask_cluster(
                num_jobs,
                condorsubmit=condorsubmit,
                condorenv=condorinit,
                logdir=self.logdir,
                memory=memory,
                runtime=runtime
            )
            self.client = dask.distributed.Client(dask_cluster)
        self.retries = retries

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def _dask_map(self, function, *iterables, key=None):
        if len(iterables) == 0:
            return

        client = self.client
        # Setting pure=False allows resubmission, with pure=True dask assumes
        # the error that happened once will always happen and skips retrying
        tasks = client.map(function, *iterables, pure=False, key=key)
        tasks_to_itemidx = dict(zip(tasks, range(len(tasks))))
        tasks = dask.distributed.as_completed(tasks)
        task_failures = {}
        for task in tasks:
            try:
                # This call should return immediately but sometimes Dask gets
                # stuck here. Unknown why. Specify timeout to circumvent.
                result = task.result(timeout=1)
            except asyncio.exceptions.TimeoutError:
                # Retry but silence the error
                tasks.add(self._dask_resubmit_failed_task(
                    function, task, tasks_to_itemidx, iterables, key))
            except Exception as e:
                logger.exception(e)
                failures = task_failures.get(task, 0)
                if self.retries is not None and failures >= self.retries:
                    raise
                logger.info(
                    f"Task failed {failures} times and will be retried")

                new_task = self._dask_resubmit_failed_task(
                    function, task, tasks_to_itemidx, iterables, key)
                task_failures[new_task] = failures + 1
            else:
                if result is None:
                    logger.error("Task returned 'None' (usually due to dask "
                                 "killing this worker).")
                    failures = task_failures.get(task, 0)
                    if self.retries is not None and failures >= self.retries:
                        raise RuntimeError(
                            "Number of retries was exceed by a task returning "
                            "'None'. This is usually due to dask killing a "
                            "worker for exceeding memory usage.")
                    logger.info(
                        f"Task failed {failures} times and will be retried")

                    new_task = self._dask_resubmit_failed_task(
                        function, task, tasks_to_itemidx, iterables, key)
                    task.cancel()
                    task_failures[new_task] = failures + 1
                else:
                    yield result
            del tasks_to_itemidx[task]
            if task in task_failures:
                del task_failures[task]

    def _dask_resubmit_failed_task(
            self, function, task, tasks_to_itemidx, iterables, key):
        idx = tasks_to_itemidx[task]
        item = (list(args)[idx] for args in iterables)
        if key is not None:
            key = key[idx] + "-retry-" + str(uuid.uuid4())
        new_task = self.client.submit(function, *item, pure=False, key=key)
        tasks_to_itemidx[new_task] = idx
        return new_task

    def process(self, function, *iterables, key=None):
        """Call function on each item in iterables, either locally or on
        HTCondor
        The parameters are handled in the same fashion as in Python's ``map``
        function.

        Parameters
        ----------
        function
            Function to be called
        *iterables
            Arguments to the function call.
        key
            Key as in ``distributed.Client.map``

        Yields
        -------
            Return values of function in the order they are finished
        """
        if self.client is None:
            for args in zip(*iterables):
                yield function(*args)
        else:
            yield from self._dask_map(function, *iterables, key=key)

    @property
    def dashboard_link(self):
        """URL to the Dask client dashboard"""
        if self.client is None:
            return None
        return self.client.dashboard_link

    @staticmethod
    def set_global_config():
        """Set the config of the local process that is needed to errorlessly
        run on HTCondor. For example ensuring the maximum number of connections
        is large enough
        """
        # Increase maximum number of connections. Dask needs ~4 per job
        nfilelimit = resource.getrlimit(resource.RLIMIT_NOFILE)[1]
        resource.setrlimit(resource.RLIMIT_NOFILE,
                           (nfilelimit, nfilelimit))

        # Increase log level of dask to hide misc messages
        logging.getLogger("distributed").setLevel(logging.WARNING)

    def close(self):
        """Close the Dask client"""
        if self.client is not None:
            self.client.close()

    @staticmethod
    def get_enumerated_dir(parentdir):
        """Get a path to a newly made directory within parentdir. """
        i = 0
        while os.path.exists(os.path.join(parentdir, str(i).zfill(3))):
            i += 1
        path = os.path.join(parentdir, str(i).zfill(3))
        os.makedirs(path)
        return path
