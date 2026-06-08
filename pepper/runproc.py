#!/usr/bin/env python3

import os
import sys
import importlib
import argparse
import logging
from datetime import datetime

import pepper
import pepper.executor
import pepper.htcondor


BUILTIN_PROCESSORS = {
    "ttbarll": ("pepper.processor_ttbarll", "Processor")
}


class _ListShortcutsAction(argparse.Action):
    def __init__(self, option_strings, dest=argparse.SUPPRESS,
                 default=argparse.SUPPRESS, help=None):
        super().__init__(option_strings=option_strings, dest=dest,
                         default=default, nargs=0, help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        shortcuts = ", ".join(BUILTIN_PROCESSORS.keys())
        print(f"Available shortcuts: {shortcuts}")
        parser.exit()


def run_processor(processor_class=None, description=None, mconly=False):
    if description is None:
        description = ("Run a processor on files given in configuration and "
                       "save the output.")
    parser = argparse.ArgumentParser(description=description)
    if processor_class is None:
        parser.add_argument(
            "processor", help="Python source code file of the processor to "
            "run or a shortcut for a processor")
        parser.add_argument(
            "--listshortcuts", action=_ListShortcutsAction, help="List the "
            "available  shortcuts for processors to use in the processor "
            "argument and exit")
    parser.add_argument("config", help="JSON configuration file")
    parser.add_argument(
        "--eventdir", help="Event destination output directory. If not "
        "specified, no events will be saved")
    parser.add_argument(
        "-o", "--output", help="Directory to save final output to. This "
        "usual are objects like cutflow and histograms. Defaults to '.'",
        default=".")
    parser.add_argument(
        "--file", nargs=2, action="append", metavar=("dataset_name", "path"),
        help="Can be specified multiple times. Ignore datasets given in "
        "config and instead process these. Can be specified multiple times.")
    parser.add_argument(
        "--dataset", action="append", help="Only process this dataset. Will "
        "also match regex expressions, though these need to be given inside "
        "quoted strings in the command line. Can be specified multiple times.")
    if not mconly:
        parser.add_argument(
            "--mc", action="store_true", help="Only process MC files. Ignored "
                                              "if --file is present")
        parser.add_argument(
            "--nomc", action="store_true", help="Skip MC files. Ignored if "
                                                "--file is present")
    parser.add_argument(
        "-c", "--condor", type=int, const=10, nargs="?", metavar="simul_jobs",
        help="Split and submit to HTCondor. By default a maximum of 10 condor "
        "simultaneous jobs are submitted. The number can be changed by "
        "supplying it to this option.")
    parser.add_argument(
        "-p", "--processes", type=int, nargs="?", help="Parallelize on the local"
        "machine on a given number of processes."
    )
    parser.add_argument(
        "-r", "--retries", default=10, type=int, help="Number of times to "
        "retry if there is exception in an HTCondor job. Default is ten.")
    parser.add_argument(
        "--chunksize", type=int, help="Number of events to "
        "process at once. A smaller value means less memory usage. Defaults "
        "to 5*10^5")
    parser.add_argument(
        "--force_chunksize", action="store_true", help="If present, makes the "
        "processor process exactly the number of events given in chunksize, "
        "unless there aren't enough events in the file. This will make "
        "reading slower.")
    parser.add_argument(
        "-d", "--debug", action="store_true", help="Enable debug messages and "
        "only process a small amount of events to make debugging fast. This "
        "changes the default of --chunksize to 10000")
    parser.add_argument(
        "-l", "--loglevel",
        choices=["critical", "error", "warning", "info", "debug"],
        help="Set log level. Overwrites what is set by --debug. Default is "
        "'warning'")
    parser.add_argument(
        "-i", "--condorinit",
        help="Shell script that will be sourced by an HTCondor job after "
        "starting. This can be used to setup environment variables, if using "
        "for example CMSSW. If not provided, the local content of the "
        "environment variable PEPPER_CONDOR_ENV will be used as path to the "
        "script instead.")
    parser.add_argument(
        "--condorsubmit",
        help="Text file containing additional parameters to put into the "
        "HTCondor job submission file that is used for condor_submit"
    )
    parser.add_argument(
        "-w", "--condorworkers", type=int, default=1,
        help="Number of workers (processes) to run in parallel inside a "
        "single HTCondor job. More workers allow faster processing but also "
        "increases memory usage inside the job. Default is 1."
    )
    parser.add_argument(
        "-m", "--memory", type=float, help="Memory in GiB that is requested "
        "per condor worker. If a worker exceeds the limit, condor might "
        " kill it. Default is 2.", default=2.0)
    parser.add_argument(
        "--runtime", type=float, help="Runtime in hours that a condor job "
        "is allowd to take. Note that jobs that exceed the runtime and are "
        "killed will still be resubmitted. Default is 3.", default=3.0)
    parser.add_argument(
        "--condorlogdir", help="Directory to store stdout and stderr logs "
        "running on HTCondor. Default is pepper_logs", default="pepper_logs")
    parser.add_argument(
        "--dasklogs",
        action="store_true",
        help="Flag to enable logging of output from dask_jobqueue. "
        "Enabling this might help with debugging condor issues, but leads "
        "to very verbose output in the console.")
    parser.add_argument(
        "--metadata", help="File to cache metadata in. This allows speeding "
        "up or skipping the preprocessing step. Default is "
        "'pepper_metadata.coffea'", default="pepper_metadata.coffea")
    parser.add_argument(
        "--statedata", help="File to write and load processing state to/from. "
        "This allows resuming the processor after an interruption that made "
        "the process quit. States produced from different configurations "
        "should not be loaded, as this can lead to bogus results. If not "
        "given, the file is put in the output directory as controlled by "
        "the --output argument.", default=None
    )
    parser.add_argument(
        "-R", "--resume", action="store_true", help="If present and the file "
        "pointed to by --statedata exists, load the file and resume from that "
        "state"
    )
    args = parser.parse_args()

    if args.condor and args.processes:
        sys.exit("Cannot give --condor and --processes at the same time")

    logger = logging.getLogger("pepper")
    logger.addHandler(logging.StreamHandler())
    if args.debug:
        logger.setLevel(logging.DEBUG)
    if args.loglevel is not None:
        if args.loglevel == "critical":
            logger.setLevel(logging.CRITICAL)
        elif args.loglevel == "error":
            logger.setLevel(logging.ERROR)
        elif args.loglevel == "warning":
            logger.setLevel(logging.WARNING)
        elif args.loglevel == "info":
            logger.setLevel(logging.INFO)
        elif args.loglevel == "debug":
            logger.setLevel(logging.DEBUG)

    if processor_class is None and args.processor in BUILTIN_PROCESSORS:
        name, class_name = BUILTIN_PROCESSORS[args.processor]
        module = importlib.import_module(name)
        Processor = getattr(module, class_name)
        logger.debug(f"Using processor from {module.__file__}")
    elif processor_class is None:
        spec = importlib.util.spec_from_file_location("pepper", args.processor)
        if spec is None:
            sys.exit(f"No such processor shortcut or file: {args.processor}")
        proc_module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(proc_module)
        try:
            Processor = proc_module.Processor
        except AttributeError:
            sys.exit("Could not find class with name 'Processor' in "
                     f"{args.processor}")
    else:
        Processor = processor_class

    config = Processor.config_class(args.config)
    store = None
    if ("store" in config
        and ("file_mode" not in config
             or "local" in config["file_mode"])):
        store = config["store"]
    if store is not None and not os.path.exists(store):
        raise pepper.config.ConfigError("store directory does not exist")

    if "file_mode" in config and "xrootd" in config["file_mode"] and (
            "X509_USER_PROXY" not in os.environ
            or not os.path.exists(os.environ["X509_USER_PROXY"])):
        raise RuntimeError(
            "xrootd usage requested but environment variable "
            "X509_USER_PROXY does not point to a proxy certificate. "
            "Please set X509_USER_PROXY and run voms-proxy-init --voms "
            "cms --out $X509_USER_PROXY")

    processor = Processor(config, args.eventdir)

    datasets = {}
    if args.file is None:
        if not config["compute_systematics"]:
            exclude = config["dataset_for_systematics"].keys()
        else:
            exclude = None
        if not mconly and args.mc and args.nomc:
            sys.exit("--mc and --nomc cannot both be present")
        datasets = config.get_datasets(
            args.dataset, exclude, "mc" if mconly or args.mc else
                                   "data" if args.nomc else "any")
    else:
        datasets = {}
        for customfile in args.file:
            if customfile[0] in datasets:
                datasets[customfile[0]].append(customfile[1])
            else:
                datasets[customfile[0]] = [customfile[1]]

    if args.file is None:
        num_files = sum(len(dsfiles) for dsfiles in datasets.values())
        num_mc_files = sum(len(datasets[dsname])
                           for dsname in config["mc_datasets"].keys()
                           if dsname in datasets)

        print(
            f"Got a total of {num_files} files of which {num_mc_files} are MC")

    if args.debug:
        print("Processing only one chunk per dataset because of --debug")
        datasets = {key: [val[0]] for key, val in datasets.items()}

    if len(datasets) == 0:
        sys.exit("No datasets found")

    # Create histdir and in case of errors, raise them now (before processing)
    os.makedirs(args.output, exist_ok=True)

    if args.statedata is None:
        statedata = os.path.join(args.output, "pepper_state.coffea")
    else:
        statedata = args.statedata
    if os.path.realpath(args.metadata) == os.path.realpath(statedata):
        print("--metadata and --statedate can not be the same")
        sys.exit(1)
    if os.path.exists(statedata) and os.path.getsize(statedata) > 0:
        if not args.resume:
            print("Found old processor state file. Please either delete "
                  f"'{statedata}' or specify --resume/-R")
            sys.exit(1)
        mtime = datetime.fromtimestamp(os.stat(statedata).st_mtime)
        print(f"Loading old processor state made on {mtime}")

    datasets = processor.preprocess(datasets)
    if args.condor is not None:
        pepper.htcondor.Cluster.set_global_config(dasklogs=args.dasklogs)
    if "exit_on_failed_jobs" in config:
        exit_on_failed_jobs = config["exit_on_failed_jobs"].lower()
    else:
        exit_on_failed_jobs = "all"
    cluster = pepper.htcondor.Cluster(
        args.processes if args.processes is not None else args.condor,
        condorsubmitfile=args.condorsubmit,
        condorinit=args.condorinit,
        retries=args.retries,
        exit_on_failed_jobs=exit_on_failed_jobs,
        mc_dsnames=list(config["mc_datasets"].keys()),
        logdir=os.path.realpath(args.condorlogdir),
        memory=str(args.memory) + " GiB",
        runtime=int(args.runtime*60*60),
        use_local=args.processes is not None
    )
    pre_executor = pepper.executor.ClusterExecutor(
        state_file_name=args.metadata,
        cluster=cluster
    )
    executor = pepper.executor.ClusterExecutor(
        state_file_name=statedata,
        cluster=cluster
    )

    try:
        pre_executor.load_state()
    except (FileNotFoundError, EOFError):
        pass
    try:
        executor.load_state()
    except (FileNotFoundError, EOFError):
        pass
    userdata = executor.state["userdata"]
    if "chunksize" in userdata:
        if (userdata["chunksize"] is not None
                and userdata["chunksize"] != args.chunksize):
            sys.exit(
                f"'{statedata}' got different chunksize: "
                f"{userdata['chunksize']}. Delete it or change --chunksize")

    maxchunks = 1 if args.debug else None
    if args.chunksize is not None:
        chunksize = args.chunksize
    elif args.debug:
        chunksize = 10000
    else:
        chunksize = 500000
    userdata["chunksize"] = args.chunksize

    if args.condor is not None:
        print(f"Dashboard available at {cluster.dashboard_link}")

    runner = pepper.executor.Runner(
        executor, pre_executor, chunksize=chunksize, maxchunks=maxchunks,
        align_clusters=not args.force_chunksize, schema=processor.schema_class,
        xrootdtimeout=pepper.misc.XROOTDTIMEOUT)
    xrootddomain = None
    if "file_mode" in config and "xrootd" in config["file_mode"]:
        xrootddomain = config["xrootddomain"]
    local_file_blacklist = None
    if "local_file_blacklist" in config \
            and config["file_mode"] == "local+xrootd":
        local_file_blacklist = config["local_file_blacklist"]
    xrootd_url_blacklist = None
    if "xrootd_url_blacklist" in config:
        xrootd_url_blacklist = config["xrootd_url_blacklist"]
    check_can_open_local = config.get("check_can_open_local_files", False)
    use_eos_redirector = config.get("use_eos_redirector", True)
    xrootd_url_priority = config.get("xrootd_url_priority", None)
    metadata = {"store_path": store, "xrootddomain": xrootddomain,
                "local_file_blacklist": local_file_blacklist,
                "url_blacklist": xrootd_url_blacklist,
                "use_eos_redirector": use_eos_redirector,
                "url_priority": xrootd_url_priority,
                "check_can_open_local": check_can_open_local,
                }
    # Give metadata for the processing step
    processor.pepperitemmetadata = metadata
    # For the preprocessing step, add metadata also to the file meta
    datasets = {k: {"files": v, "metadata": metadata}
                for k, v in datasets.items()}
    output = runner(datasets, "Events", processor)
    processor.save_output(output, args.output)
    if args.eventdir is not None:
        print("If there were errors, please make sure to run "
              "delete_duplicate_outputs.py for the per-event output in "
              f"{args.eventdir}")

    return output


if __name__ == "__main__":
    run_processor()
