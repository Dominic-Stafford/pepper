Conceptual Overview
===================

Pepper is a multi-purpose framework designed for performing event selection and data reduction of CMS nanoAOD datasets. 
It builds upon the coffea and awkward-array libraries and thus serves as a python-first alternative to traditional C++ based analysis frameworks.

The Coffea Paradigm
--------------------

To understand how Pepper operates, it is instructive to first understand coffea and the idea of columnar data analysis. 
In high-energy physics (HEP), datasets are often structured as collections of events, where each event contains multiple particles and associated attributes.
In an event-based analysis framework, one typically processes each event sequentially, applying selection criteria and computing derived quantities on an event-by-event basis.
This approach can be inefficient, especially for large datasets, as it may involve significant overhead in terms of data loading and function calls. 
Coffea, on the other hand, adopts a columnar data analysis paradigm.
Instead of processing events one at a time, coffea operates on entire columns of data simultaneously. 
In fact, coffea leverages the Awkward Array library to handle variable-length arrays efficiently. Under the hood, 
Awkward Array stores events in memory using a columnar format, where each nanoAOD attribute of particles 
is stored in a separate contiguous array. This means that selection criteria and computations are applied to entire 
arrays of particle attributes at once, leveraging the power of vectorized operations provided by libraries like NumPy and Awkward Array.
This columnar approach can lead to significant performance improvements, as it reduces the overhead associated with event-by-event processing and allows for more efficient use of modern CPU architectures.

In Coffea, analyses are structured around a central ``Processor`` class that defines the analysis workflow with the 
analysis implemented in a central ``process`` method that operates on NanoEvents objects. Delivery of events to the processor is 
handled by the ``Runner`` class, which supports both local and distributed execution. This setup allows for scaling analyses 
from a single machine to large computing clusters with minimal overhead.

Pepper's Architecture
---------------------

Pepper extends the coffea framework in numerous ways to provide a more user-friendly and feature-rich experience for HEP analysts. 
First, a ``Selector`` object handles all the bookkeeping associated with selection cuts, derived quantity definitions, categorizations, corrections
and systematic uncertainty weights. Any update of the selector object by one of its setter or adder methods automatically triggers
filling of outputs defined in the config-file by the ``OutputFiller`` class which also keeps track of the socalled "cutflows"
- tables that summarize how many events pass each selection cut.

Second, the abstract Coffea ``Processor`` class is subclassed by Pepper's ``Processor`` which handles as the technical interface 
between the ``Selector`` and the Coffea framework. In turn, the ``Processor`` is subclassed by the ``ProcessorBasicPhysics`` class which implements
standard HEP analysis features such as POG recommended cuts, object cleaning, applying corrections, trigger selection, and MET filters.
The user can then implement their analysis by subclassing ``ProcessorBasicPhysics`` and need only to implement the `process_selection` method
defining the analysis-specific selection. The figure below illustrates the class hierarchy of Pepper's processor classes and the relationship 
to Coffea's processor class.

.. image:: /img/banana.png

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The ``process_selection`` Method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The signature of the ``process_selection`` method is as follows:

.. code-block:: python3

    def process_selection(self,
                          selector: Selector,
                          dsname: str,
                          is_mc: bool,
                          filler: OutputFiller) -> None:
        ...

This method is called for each dataset (``dsname``) being processed. The ``is_mc`` boolean indicates whether the dataset is Monte Carlo simulation or real data
since different selections or corrections may apply. Through the ``selector`` object, the user has access to all event data and 
can apply selection cuts, define derived quantities. Typically, the user will not need to interact directly with the ``filler`` object,
as the ``selector`` handles output filling automatically when cuts or derived quantities are defined.

A simple example implementation of the ``process_selection`` method which is found in the ``examples/example_processor.py`` example script is as follows:

.. code-block:: python3

    def process_selection(self, selector, dsname, is_mc, filler):
        # Implement the selection steps: add cuts, define objects and/or
        # compute event weights

        # Add a cut only allowing events according to the golden JSON
        # The good_lumimask method is specified in pepper.ProcessorBasicPhysics
        # It also requires a lumimask to be specified in config
        if not is_mc:
            selector.add_cut("Lumi", partial(
                self.good_lumimask, is_mc, dsname))

        # Only allow events that pass triggers specified in config
        # This also takes into account a trigger order to avoid triggering
        # the same event if it's in two different data datasets.
        pos_triggers, neg_triggers = pepper.misc.get_trigger_paths_for(
            dsname, is_mc, self.config["dataset_trigger_map"],
            self.config["dataset_trigger_order"])
        selector.add_cut("Trigger", partial(
            self.passing_trigger, pos_triggers, neg_triggers))

        # Pick electrons satisfying our criterias
        selector.set_column("Electron", self.pick_electrons)
        # Also pick muons
        selector.set_column("Muon", self.pick_muons)

        # Only accept events that have to leptons
        selector.add_cut("Exactly2Leptons", self.lepton_pair)

        # Only accept events that have oppositely changed leptons
        selector.add_cut("OCLeptons", self.opposite_sign_lepton_pair)

        # Things that could be done next: Scaling the leptons according to POG recipes,
        # adding cuts on the jets, MET or only allowing events that have a certain m_ll.
        # A full implementation can be found in processor_ttbarll.py.

This example shows how we apply the golden JSON selection for data, require events to pass certain triggers,
select electrons and muons based on criteria defined in the ``config.json``, and apply cuts to ensure events have exactly two oppositely charged leptons.
Notice, the interplay between defining object selections using ``set_column`` and applying event-level cuts 
- on newly defined columns - using ``add_cut``. Furthermore, since we inherit from the ``ProcessorBasicPhysics`` class,
we get access to several utility methods such as ``good_lumimask`` and ``pick_electrons`` that implement common HEP analysis tasks.
Most cuts can be configured through the config file which is passed to the processor upon initialization and the config file can be 
further customized to suit specific user needs as described in detail on [this] page.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Defining Columns and Applying Cuts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Running Your Processor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Use ``pepper.runproc``. Refer to this other page for a more detailed guide on how to run your analysis.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Defining Histograms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
