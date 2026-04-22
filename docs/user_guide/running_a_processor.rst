.. _running-your-processor:

Running Your Processor
======================

Once you have defined a physics processor for your analysis, pepper needs to know in what context to run it. 
This includes where to run (locally or on a cluster), what data to use (which config file), and where to output the results set. 
The command-line program ``pepper.runproc`` is used to take care of this. Overall, the command looks like this:

.. code-block:: bash

   python -m pepper.runproc <custom_processor.py> --config <path_to_config_file> --output <output_directory>


When executed this way, all the datasets defined in the config file will be passed to the processor for local processing 
and aggregated histograms and cutflows will be saved to the output directory. Typically, this does not reflect how we
would use the processor. In a normal HEP workflow, we would initially run the processor on a subset of the data to test it. 
When we are happy with the results (typically when we see no errors), we would then run it on the full dataset - typically on a cluster.
The ``pepper.runproc`` script supports both these use-cases. For processing only the first chunk of each dataset, pass the debug flag ``--debug``.
If you want to submit to HTCondor, just add the ``--condor WORKERS`` argument and specify the number of workers.

.. note::

    It is also possible to submit to HTCondor when running in debug mode. This feature is typically used when testing whether the systematics calculation works as expected.

To control which environment is employed on the HTCondor node, the parameter ``--condorinit`` can be used. 
``--condorinit`` should point to a Shell script that can be sourced setting up the environment. 
If ``--condorinit`` is not present, Pepper will instead use the script that is pointed at by the local environment 
variable ``PEPPER_CONDOR_ENV``. If this is also not set, the jobs will be run in the default environment of your HTCondor system.


.. note::

    The terminal running the Pepper processor needs to stay alive until the processing is complete.
    This means you should avoid closing the terminal or interrupting the process until you see a completion message.
    If the process is interrupted, you may need to restart it from the beginning. On NAF, this can be done using the
    ``tmux``, ``nohup`` or ``screen`` commands to run the process in the background while it has to be ``tmux`` on LXPLUS. 
    If you haven't already, ``tmux`` must be enabled on LXPLUS using ``systemctl --user start tmux.service && tmux a``.

A directory with logs from the jobs will be present under ``pepper_logs``. Directories inside ``pepper_logs`` are numbered, 
the highest number is the one of the latest run. The log level of the logs inside is controller via the ``--loglevel`` option. 
Set it to ``debug`` to get full logging.

If you have access to a dedicated machine, it is also possible to run pepper locally on multiple CPU cores by giving the ``--processes WORKERS`` argument.
Please do not abuse this feature on the login nodes of clusters!

Further options can be seen by running ``python -m pepper.runproc -h``.

When Running a Processor Fails
-------------------------------

If the processor fails partway through, it does not need to be restarted from
scratch. Resume processing by running the same command again with the
``--resume`` flag:

.. code-block:: bash

   python -m pepper.runproc <custom_processor.py> --config <path_to_config_file> --output <output_directory> --resume

Once the processor has completed successfully, run the cleanup script to remove
any duplicate outputs that may have been produced across the original and
resumed runs:

.. code-block:: bash

   python scripts/delete_duplicate_outputs.py