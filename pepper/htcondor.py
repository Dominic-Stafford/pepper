import os
import logging
import resource
import dask
import dask.distributed
import dask_jobqueue
import shlex
import traceback
from collections import defaultdict

import pepper.config

logger = logging.getLogger(__name__)


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


def get_dask_cluster(num_jobs, runtime=3*60*60, memory="2 GiB", disk="3 GiB",
                     cores=1, *, condorsubmit=None, condorenv=None,
                     logdir=None):
    """Get a Dask Jobqueue HTCondor cluster for a host

    Parameters
    ---------
    num_jobs
        Number of jobs to run in parallel
    runtime
        Requested runtime in seconds. If None, do not request a runtime
    memory
        Request memory. String with a unit like "GiB" or an int.
    disk
        Request disk space. String with a unit like "GiB" or an int.
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
        "+RequestRuntime": str(int(runtime)),
        "+MaxRuntime": str(int(runtime))
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
        # Set memory limit to 0 as explained above
        memory=0,
        disk=disk,
        log_directory=logdir,
        job_extra_directives=job_extra_directives,
        job_script_prologue=["source " + shlex.quote(condorenv)],
        # Set port parameters to 0 to use random ports
        scheduler_options={"dashboard_address": ":0"}
    )
    site = get_site()
    config.update(site_config.get(site, {}))
    cluster = dask_jobqueue.htcondor.HTCondorCluster(**config)
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


class PepperSchedulerPlugin(
        dask.distributed.diagnostics.plugin.SchedulerPlugin):
    """A dask scheduler plugin which prints certain events - adding workers,
    removing workers, and task state transitions - to the default log for
    better diagnostics.
    See https://distributed.dask.org/en/latest/plugins.html for the plugin API.
    """

    def __init__(self):
        self.current_workers = []

    def add_worker(self, scheduler, worker, **kwargs):
        """Called when the dask scheduler successfully connects to a worker.
        """

        if worker in scheduler.workers:
            worker_name = scheduler.workers[worker].name
        else:
            worker_name = "unknown"
        logger.info(f"Connected to worker at adress {worker} "
                    f"with name {worker_name}")
        self.current_workers.append(worker)

    def remove_worker(self, scheduler, worker, **kwargs):
        """Called when a worker is removed from the scheduler for any reason
        (including both dead workers and regular shutdowns).
        """

        # Only log if the worker is expected to still run
        if len(scheduler.tasks) > 0 and worker in self.current_workers:
            logger.warn(f"HTCondor worker process at address {worker} died. "
                        "Resubmitting automatically.")
            self.current_workers.remove(worker)

    def transition(self, key, start, finish, *args, **kwargs):
        """Called when a task changes state from 'start' to 'finish'.
        There are quite a few different states, we only log the most important
        """

        if finish == "processing":
            logger.info(f"Started processing task {key}")
        elif start == "processing" and finish == "memory":
            if "startstops" in kwargs and len(kwargs["startstops"]) > 0:
                startstops = kwargs["startstops"][-1]
                time = startstops["stop"] - startstops["start"]
                logger.info(f"Finished processing task {key} "
                            f"in {time:.1f} seconds")
            else:
                logger.info(f"Finished processing task {key}")
        elif start == "processing" and finish != "memory":
            logger.info(f"Task {key} did not finish processing "
                        f"(status '{finish}')")

    def valid_workers_downscaling(self, scheduler, workers):
        """Called when Dask decides that some workers should be shut down
        since there are no more tasks for them to process
        This can in principle veto the shutdown, but we do not need that
        """

        for worker in workers:
            address = worker.address
            if address in self.current_workers:
                logger.debug(f"Shutting down the worker at address {address}")
                self.current_workers.remove(address)
        # Return all workers, i.e. allow the shutdown of all of them
        return workers


class Cluster:
    """Run tasks either locally or on HTCondor.

    Internally, this class may call ``get_dask_cluster`` in order to submit
    jobs to HTCondor.
    """
    def __init__(
            self, num_jobs, condorsubmit=None, condorinit=None,
            logdir="pepper_logs", retries=None, exit_on_failed_jobs="all",
            mc_dsnames=[], condorsubmitfile=None, memory="2 GiB",
            runtime=3*60*60):
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
        exit_on_failed_jobs
            Whether task which fail retries times should cause an exception.
            Can be "All" if all tasks should cause a failure, "Data" if only
            data tasks should cause an exception, or "None" to never raise the
            exception
        mc_dsnames
            A list of the MC dataset names. Needed for exit_on_failed_jobs
            == "data"
        condorsubmitfile
            Path to a file containing additional content to add to the
            HTCondor submit file
        memory
            Request memory. String with a unit like "GiB" or an int.
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
            self.client.register_plugin(PepperSchedulerPlugin())
        self.retries = retries
        if exit_on_failed_jobs not in ["all", "data", "none"]:
            raise pepper.config.ConfigError(
                "'exit_on_failed_jobs' must be one of 'all', 'data' or 'none'")
        self.exit_on_failed_jobs = exit_on_failed_jobs
        self.mc_dsnames = mc_dsnames

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def _dask_map(self, function, *iterables, key=None):
        if len(iterables) == 0:
            return

        client = self.client
        tasks = client.map(function, *iterables, pure=True, key=key)
        # Get an iterator that yields tasks in the order they complete
        tasks_iterator = dask.distributed.as_completed(tasks)
        # Dictionary to store the number of retries per task
        task_failures = defaultdict(int)
        # Store completed tasks to check for completeness at the end
        tasks_completed = []

        # Iterate over all still pending tasks
        for task in tasks_iterator:
            if not task.done():
                logger.critical(f"Task {task.key} was yielded from the "
                                "iterator even though it was not done. "
                                "This Should Not Happen (TM). "
                                "Please contact the Pepper developers.")
            if task.status == "finished":
                # Tasks that are finished successfully
                result = task.result()
                # Sanity check: the task should not be completed twice
                if task.key in tasks_completed:
                    logger.critical(f"Task {task.key} was evaluated twice! "
                                    "This Should Not Happen (TM). "
                                    "Please contact the Pepper developers. "
                                    "The duplicate result will be ignored.")
                # Check whether the result is actually there
                elif result is None:
                    logger.critical(f"Result for task {task.key} is None. "
                                    "This Should Not Happen (TM). "
                                    "Please contact the Pepper developers. "
                                    "The task will be resubmitted.")
                    task_failures[task.key] += 1
                    task.retry()
                    tasks_iterator.add(task)
                else:
                    logger.debug(f"Got result for task {task.key}")
                    tasks_completed.append(task.key)
                    # We cancel the task after we got the result so it doesnt
                    # re-run if its worker dies
                    task.cancel()
                    # Only if we got here, we further pass the result to the
                    # accumulator
                    yield result
            else:
                # Task did not finish successfully. Try to get information why
                exc = task.exception()
                tb = task.traceback()
                task_retries = task_failures[task.key]
                logger.error(f"Task failed with status '{task.status}' "
                             f"for '{task.key}' (retry {task_retries}).")
                if exc is not None:
                    logger.error(f"The type of the exception is "
                                 f"'{type(exc).__name__}'.")
                    logger.error("The exception message is: ")
                    logger.error(str(exc))
                else:
                    logger.error("The type of the exception is not available.")
                if tb is not None:
                    logger.error("Stacktrace of the exception: ")
                    traceback.print_tb(tb)
                else:
                    logger.error("No stack trace is available.")
                is_mc = (task.key.split("/")[0] in self.mc_dsnames)
                if self.retries is None or task_retries < self.retries:
                    # Retry the task using the task.retry() method
                    # We do it this way as opposed to the automatic retry
                    # functionality from Dask to log the retries, and keep
                    # more control over when a task should be retried
                    logger.error("Retrying the task automatically.")
                    task_failures[task.key] += 1
                    task.retry()
                    # Re-add the task to the iterator to wait for it again
                    tasks_iterator.add(task)
                elif (self.exit_on_failed_jobs == "all" or (
                        self.exit_on_failed_jobs == "data" and not is_mc)):
                    logger.error("Maximum number of retries reached. "
                                 "Aborting.")
                    raise exc
                else:
                    logger.error(f"Maximum number of retries reached for "
                                 f"task {task.key}. Will continue running, "
                                 f"but not retry this task.")

        logger.debug("All tasks processed. Checking for completeness...")
        incomplete_tasks = [task.key for task in tasks
                            if task.key not in tasks_completed]
        if len(incomplete_tasks) > 0:
            logger.critical(f"WARNING: {len(incomplete_tasks)} tasks that "
                            "were not processed properly! "
                            "The following tasks where not processed:")
            logger.critical('\n'.join(incomplete_tasks))
            logger.critical("This may cause issues if observed data is "
                            "skipped, or if using pre-computed lumifactors")
        else:
            logger.debug("All tasks are complete.")

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
    def set_global_config(dasklogs=False):
        """Set the config of the local process that is needed to errorlessly
        run on HTCondor. For example ensuring the maximum number of connections
        is large enough
        """
        # Increase maximum number of connections. Dask needs ~4 per job
        nfilelimit = resource.getrlimit(resource.RLIMIT_NOFILE)[1]
        resource.setrlimit(resource.RLIMIT_NOFILE,
                           (nfilelimit, nfilelimit))

        if dasklogs:
            # Adapt level for distributed logger
            logging.getLogger("distributed").setLevel(logging.DEBUG)
            # Enable and format dask_jobqueue logger for console output
            jobqueue_logger = logging.getLogger("dask_jobqueue.core")
            stream_handler = logging.StreamHandler()
            stream_handler.setFormatter(
                logging.Formatter(fmt="[%(name)s(%(levelname)s)] %(message)s"))
            jobqueue_logger.addHandler(stream_handler)
            jobqueue_logger.setLevel(logging.DEBUG)
        else:
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
