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

For a full run, submit to HTCondor by adding the ``--condor WORKERS`` argument; see :ref:`htcondor` for the complete
setup, including environment scripts, resource requests, logs, and recovery from failures.

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