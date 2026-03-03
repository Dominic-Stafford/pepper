Getting Started
===============

First, install the ``pepper`` package by following the instructions in the
:ref:`installation` section. Afterwards, you should have a working Pepper
installation with a local virtual environment set up.

In Pepper, an analysis is implemented as a ``Processor`` class. A short example
of such a processor with explanatory comments can be found in
:repo:`example/example_processor.py`. It demonstrates a basic dileptonic
selection in a :math:`t\bar{t}` event, including golden JSON and trigger requirements, electron and muon
object selection, and cuts requiring exactly two oppositely charged leptons.
The example directory also contains a corresponding :repo:`example/example_config.json`
which configures the object selections and datasets used by the processor.

This processor can be run by executing the following command from inside the
example directory:

.. code-block:: bash

   python -m pepper.runproc example_processor.py example_config.hjson

When developing or testing, pass the ``--debug`` flag to process only the first
chunk of each dataset rather than running over the full input:

.. code-block:: bash

   python -m pepper.runproc example_processor.py example_config.hjson --debug

A full list of available command line options can be seen by running:

.. code-block:: bash

   python -m pepper.runproc -h

After a successful run, the output directory will contain aggregated histograms
and a cutflow table summarising how many events passed each selection step.

For the CMS Data Analysis School (DAS) in Hamburg in 2025, one of the long
exercises used Pepper to perform a measurement of the top quark pair production
cross section using :math:`{\sim}1\;\text{fb}^{-1}` of data. This tutorial is
`available online <https://cmsdas-2025-hamburg-longex.docs.cern.ch/ex_topxs/>`__.

Feel free to ask any questions in our :mattermost:`Pepper Mattermost <pepper-users>` chat!

Next Steps
----------

Once you have the example running, the following resources are good places to
continue:

- :ref:`conceptual-overview` - explains the architecture of Pepper's
  ``Processor``, ``Selector``, and ``OutputFiller`` classes in depth.
- :repo:`pepper/processor_ttbarll.py` - a full reference implementation of a
  dileptonic :math:`t\bar{t}` analysis, including systematics and kinematic
  reconstruction.
- :ref:`configuration-reference` - a complete description of all available
  configuration options.
- :ref:`running-your-processor` - covers HTCondor submission and how to resume a failed run.