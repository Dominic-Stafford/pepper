.. _htcondor:

Running on HTCondor
===================

For anything larger than a small debug or test run, it is recommended that Pepper distributes work across an HTCondor batch system. 
This page covers the practical end of that: how to invoke condor submission, how to set up the environment that jobs run in, 
where to find logs, what to do when jobs fail, and the X509/VOMS proxy setup needed for remote-data access from inside jobs.

The narrative material on ``pepper.runproc`` itself is in :ref:`running-your-processor`. This page assumes you can already run 
the processor locally and want to scale up.


Submitting jobs
---------------

To submit, add ``--condor`` (alias ``-c``) to your ``pepper.runproc`` command. 
The optional integer argument sets the maximum number of simultaneous Condor jobs:

.. code-block:: bash

   python -m pepper.runproc my_processor.py my_config.json \
       --condor 200 --output prod_v1

Without an integer, ``--condor`` defaults to 10. The number is an upper bound on simultaneous jobs, 
not on total jobs -- Pepper splits the work into many small tasks and submits them as the cluster has capacity. 
Higher values give faster wall-clock time at the cost of more peak load on the shared system.

The Python process controlling the run stays alive until all jobs finish, so use ``tmux`` 
(or ``nohup`` / ``screen`` on systems that allow them) to survive a closed terminal:

.. code-block:: bash

   tmux new -s pepper
   # ... source environment, run the command above ...
   # detach with Ctrl-B then D; reattach with `tmux a -t pepper`

On lxplus, ``tmux`` is gated behind a per-user systemd unit and must be enabled once per account:

.. code-block:: bash

   systemctl --user start tmux.service && tmux a

If you have a dedicated multi-core machine and want to parallelise locally instead of submitting to Condor, use ``--processes WORKERS``.
**Do not do this on a shared login node** -- it will consume all available cores, and you will earn the immediate attention of the site administrators.


Environment setup inside jobs
-----------------------------

A Condor worker starts in a clean shell with none of your local modules, virtual environments, or grid tools loaded. Pepper provides two hooks to set this up.

``--condorinit`` (alias ``-i``)
   Path to a shell script that is sourced inside each worker before ``pepper.runproc`` starts processing. 
   This is where you load CMSSW, activate a virtual environment, source the LCG view, set ``X509_USER_PROXY``, and so on.

``PEPPER_CONDOR_ENV``
   Environment variable that does the same thing, used when ``--condorinit`` is not given. 
   This is automatically in :repo:`example/environment.sh` and is a convenient way to set the 
   init script without having to pass it on the command line every time.
   In case you are not using the recommended environment setup, set it in your shell rc file once and forget about it:

   .. code-block:: bash

      export PEPPER_CONDOR_ENV=/afs/cern.ch/user/<x>/<username>/pepper_env.sh

If neither is set, jobs run in whatever default environment the Condor system provides which will not have the necessary dependencies loaded. 
After all, Pepper itself, its dependencies, and the CMS grid environment all have to be made available somehow.

The :repo:`example environment script <example/environment.sh>` shows a working setup for DESY NAF and lxplus that sources the grid environment, 
loads an LCG view, activates a virtual environment, and sets the proxy path. 
It is a useful starting point even if your analysis ultimately needs something different.

.. tip::

   The init script runs on every worker, every job. Keep it fast.
   Slow operations -- compiling code, fetching from the network -- multiply across hundreds of jobs and can dominate the run.


Resource requests
-----------------

Three options control what each job asks Condor for:

``--memory N`` (alias ``-m``)
   Memory per worker in GiB. Default ``2.0``. Workers that exceed this may be killed by Condor.

``--runtime H``
   Maximum runtime per job in hours. Default ``3.0``. 
   Jobs that exceed the runtime are killed by Condor but will be automatically resubmitted by Pepper as long as retries remain.

``--condorworkers N`` (alias ``-w``)
   Number of dask workers (parallel processes) inside a single Condor job. 
   Default ``1``. Higher values run more processing in parallel per slot but raise the per-job memory request correspondingly. 
   Useful when the cluster's per-slot memory limit is high but the per-user job count limit is low.

For finer-grained control, ``--condorsubmit`` accepts a path to a file containing additional lines that are appended verbatim to the HTCondor submit description. 
Use this for site-specific knobs that Pepper does not expose directly -- requesting specific accounting groups, GPU resources, special queues, and so on.


Logs and the ``pepper_logs`` directory
--------------------------------------

By default, every run creates a new subdirectory under ``pepper_logs/`` containing the stdout, stderr, and Condor logs for that run's jobs. 
The subdirectory name is a zero-padded sequence number: the very first run is ``pepper_logs/000``, the next ``pepper_logs/001``, and so on. 
The highest-numbered subdirectory is always the most recent.

.. code-block:: text

   pepper_logs/
   |-- 000/
   |   |-- worker-0.out
   |   |-- worker-0.err
   |   |-- worker-0.log
   |   |-- ...
   |-- 001/
   |-- 002/
   |-- ...

Pepper never deletes old log directories. They accumulate and eventually become large; cleaning them up periodically is your responsibility.

The verbosity inside the logs is controlled by ``--loglevel`` 
(or ``--debug``, which sets it to ``debug`` and enables a few additional test paths). 
For diagnosing job-side failures, set ``--loglevel debug``; the resulting logs are large but include every event chunk boundary and exception traceback.

To put logs somewhere other than ``./pepper_logs``, use ``--condorlogdir``. 
This is mainly useful on lxplus, where AFS quotas make it convenient to point logs at ``/eos/user/...``.

.. note::

   If your ``--condorlogdir`` is on EOS at lxplus (``/eos/user/...``), Pepper uses :mod:`dask_lxplus`'s ``CernCluster`` instead of the standard ``HTCondorCluster``. 
   Install ``dask-lxplus`` (``pip install dask-lxplus``) before using EOS log paths there.

For diagnosing dask scheduling issues rather than analysis failures, the ``--dasklogs`` flag enables dask-jobqueue's own logging. 
This is verbose and useful primarily when jobs are failing to start at all.


Resuming after a failure
------------------------

If a run fails partway through -- because of a cluster issue, an exhausted retry budget, or because you cancelled it --
**do not start over**. Re-run the same ``pepper.runproc`` command with the ``--resume`` (alias ``-R``) flag:

.. code-block:: bash

   python -m pepper.runproc my_processor.py my_config.json \
       --condor 200 --output prod_v1 --resume

Pepper persists processing state between runs in a state file 
(by default inside the output directory, controllable with ``--statedata``). 
With ``--resume``, only chunks that did not complete successfully are re-submitted; 
completed work is loaded from the existing state and skipped.

.. caution::

   The state file is tied to a specific configuration. Resuming with a different config, different processor, or different dataset
   is undefined behaviour and will silently produce bogus results. 
   If you have changed anything other than the resume flag itself, start a fresh output directory.

Per-job retries are controlled independently by ``--retries N`` (default 10). 
This handles individual flaky workers without exiting the run. 
``--retries`` is for transient failures (e.g. node crashes); ``--resume`` is for restarting an entire interrupted run.


X509 proxy for remote data
--------------------------

If your config reads NanoAOD files over XRootD (the typical case when running on remote Tier 2 data), every Condor worker needs a VOMS proxy. 
The full XRootD setup is documented in :doc:`remote_data`; the Condor-specific points are:

- ``X509_USER_PROXY`` must point at a path that is **readable from inside the Condor worker**. 
  ``/tmp/`` is per-machine and not shared, so the default proxy location does not work. 
  A path under your home directory (e.g. ``~/.globus/x509up``) does.

- The init script (``--condorinit`` or ``PEPPER_CONDOR_ENV``) must export ``X509_USER_PROXY`` so that workers see the same value the submitter does. 
  The :repo:`example environment script <example/environment.sh>` does this.

- The proxy itself must be created once before submission:

  .. code-block:: bash

     voms-proxy-init --voms cms --out $X509_USER_PROXY

  Proxies expire (by default after 24 hours). 
  Re-run the command when you see ``sslv3 alert certificate expired`` errors or when long runs start failing partway through. 
  For runs that will outlive a single proxy, request a long-lived proxy with ``--valid 192:00`` (8 days, the maximum on most VOMS servers).


Site-specific behaviour
-----------------------

Pepper detects the site automatically from the hostname and adjusts its Condor submission accordingly. 
The two sites with built-in handling are:

**lxplus** (CERN)
   Sets the dask scheduler port to ``8786`` and constrains worker ports to ``10000-10100`` to fit through the lxplus firewall.
   Switches to ``dask_lxplus.CernCluster`` when logs are on EOS.
   Disables Condor stdout/stderr streaming (which conflicts with the batch system's own output handling).

**naf** (DESY)
   No special Condor configuration beyond the defaults, but recognised so site-aware scripts 
   (the example ``environment.sh``, future site-specific behaviours) can branch on it.

Other sites fall through to the default behaviour, which works on most generic HTCondor installations. 
If your site needs different defaults to function at all, see ``pepper.htcondor.get_dask_cluster`` in the API reference 
-- the per-site configuration lives in a small dictionary near the top of that function and is the right place to add support for a new site.


Failure handling policy
-----------------------

When the retry budget for a job is exhausted, Pepper's behaviour depends on whether the failed job was processing data or MC. 
The ``exit_on_failed_jobs`` parameter (passed through from the configuration, not directly from the command line) takes three values:

``"all"`` (default)
   Any unrecoverable job failure raises an exception and exits the run. 
   Safest default; you find out immediately when something is wrong.

``"data"``
   Only failures in *data* jobs exit the run. MC failures are tolerated. 
   Useful for long MC productions where individual MC files may have known issues that should not block the data analysis.

``"none"``
   Never exit on job failures. The run continues with whatever completed, and the output will be missing chunks. 
   Use only when you have a clear reason to allow only partial results.

In all three modes, the offending jobs' logs remain in ``pepper_logs/`` for inspection.


See also
--------

- :ref:`running-your-processor` -- the general processor invocation, output structure, and ``--resume`` semantics.
- :doc:`remote_data` -- full XRootD setup, including the priority list of fallback sites.
- :ref:`bundled-scripts-example` -- where Condor submission fits into a complete analysis workflow.
- :repo:`example/environment.sh` -- a working init script for DESY NAF and lxplus, useful as a template.