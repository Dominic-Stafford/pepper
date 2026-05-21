.. _conceptual-overview:

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

.. mermaid::

   classDiagram
       direction TB

       class `coffea.ProcessorABC` {
           <<abstract>>
       }
       class `pepper.Processor`
       class `pepper.ProcessorBasicPhysics`
       class UserProcessor {
           +process_selection()
       }
       class `pepper.Selector`
       class `pepper.OutputFiller`

       `coffea.ProcessorABC` <|-- `pepper.Processor`
       `pepper.Processor` <|-- `pepper.ProcessorBasicPhysics`
       `pepper.ProcessorBasicPhysics` <|-- UserProcessor
       `pepper.Processor` *-- `pepper.Selector`
       `pepper.Selector` ..> `pepper.OutputFiller` : update triggers fill

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The ``process_selection`` Method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The signature of the ``process_selection`` method is as follows:

.. code-block:: python

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

.. code-block:: python

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

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Defining Columns and Applying Cuts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The two primary operations within ``process_selection`` are defining columns on the
selector and adding cuts. These correspond to object selection and event-level selection
respectively.

""""""""""""""""""""
Defining columns
""""""""""""""""""""

Object selections are imposed by setting columns on the selector using ``set_column``.
For example, to apply standard electron selections::

    selector.set_column("Electron", self.pick_electrons)

This tells the selector to set - in this case overwriting - the column named
``"Electron"`` using the function ``pick_electrons``, which is called with all data
currently in the selector. Users can define this function themselves for non-standard
object definitions, but for standard POG-recommended selections the function is already
implemented in ``ProcessorBasicPhysics``. In that case, the user only needs to specify
the desired cuts in the configuration file:

.. code-block:: json

    "ele_cut_transreg": true,
    "ele_eta_min": -2.4,
    "ele_eta_max": 2.4,
    "good_ele_id": "mva:Iso90",
    "good_ele_pt_min": 20.0,

""""""""""""""""""""
Adding cuts
""""""""""""""""""""

Event-level cuts are applied using ``add_cut``::

    selector.add_cut("AtLeast2Leps", partial(self.lepton_pair, is_mc))

Here ``"AtLeast2Leps"`` is the name of the cut, implemented by the ``lepton_pair``
function. The function is also passed a flag indicating whether the current dataset is
Monte Carlo simulation or observed data, since different selections or corrections may
apply to each.

A cut function can return either a boolean array, where ``False`` entries are removed,
or a numeric array, where zero entries are discarded and all surviving events are scaled
by the corresponding value. Pepper uses this convention to encapsulate scale factor
corrections together with the selection step they correct for. For example:

.. code-block:: python

    def lepton_pair(self, is_mc, data):
        """Select events that contain at least two leptons."""
        accept = np.asarray(ak.num(data["Lepton"]) >= 2)
        if is_mc:
            weight, systematics = self.compute_lepton_sf(data[accept])
            accept = accept.astype(float)
            accept[accept.astype(bool)] *= np.asarray(weight)
            return accept, systematics
        else:
            return accept

For MC events, the systematic weight variations corresponding to the lepton scale factors
are returned as a dictionary alongside the selection array. If systematics computation is
disabled in the config, this dictionary will be empty.

""""""""""""""""""""
Categorisations
""""""""""""""""""""

Analyses often define different event categories with slightly different selections. Once
boolean arrays have been computed for each category - for instance the lepton-flavour
channels in a :math:`t\bar{t}` analysis - a categorisation can be registered as:

.. code-block:: python

    selector.set_cat("channel", {"is_ee", "is_em", "is_mm"})

where ``"channel"`` is the name of the categorisation. All histograms produced will by
default be split according to this categorisation. A special case is a categorisation named
``"dataset"``: if set, it is used in place of the dataset name when filling histograms.
This is particularly useful for unfolding measurements that need to split signals into
different particle-level bins.

""""""""""""""""""""""""""""""""""""""""
Jet energy scale systematics
""""""""""""""""""""""""""""""""""""""""

One set of systematic variations that cannot be handled as simple event weights are jet
energy scale variations, since all selection steps involving jets must be re-performed for
each variation. It is recommended to implement these in a separate function
``process_selection_jet_part``, which can then be called once per variation. Pepper
provides utilities for rescaling jets according to each variation. Other sample-based
systematics need only be declared in the config file, and Pepper will record them in
histograms as a systematic entry on the corresponding dataset.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Defining Histograms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Histograms are defined declaratively in the HJSON configuration file rather than constructed in processor code. 
Each entry in the ``hists`` object names a histogram, lists its axes, 
and points at the per-event data to fill it from using Pepper's ``DataPicker`` syntax. 
By default Pepper produces one histogram per cut, so the cutflow and histogram outputs stay in sync without extra user effort.

See :ref:`histograms` for the full lifecycle 
-- defining histograms, common DataPicker patterns, the output file layout, and how to read histograms back into Python.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Deferring Cut Application
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, every call to ``add_cut`` immediately discards events that fail the cut from
the selector's internal data array. This behaviour can be suspended by setting:

.. code-block:: python

    selector.applying_cuts = False

While ``applying_cuts`` is ``False``, calls to ``add_cut`` store the cut in an internal
queue (``unapplied_cuts``) rather than removing any rows. The full event array is
preserved, and the combined effect of all queued cuts is accessible through two
properties:

- ``selector.final`` - the data array masked to only the events that pass all queued cuts.
- ``selector.final_systematics`` - the corresponding systematics, masked in the same way.
- ``selector.num_final`` - the count of events passing all queued cuts.

Callbacks triggered by ``add_cut`` still fire during this mode, but they receive
``selector.final`` and ``selector.final_systematics`` rather than the full array, so they
see only the would-be-selected events.

""""""""""""""""""""""""""""""""""""""""
Re-enabling cut application
""""""""""""""""""""""""""""""""""""""""

When ``applying_cuts`` is set back to ``True``, any accumulated unapplied cuts are
flushed immidiately - all events failing the queued cuts are discarded at once:

.. code-block:: python

    selector.applying_cuts = True   # triggers apply_all_cuts() if cuts are pending

""""""""""""""""""""""""""""""""""""""""
Defining columns under deferred cuts
""""""""""""""""""""""""""""""""""""""""

The ``set_column`` method accepts an ``all_cuts`` keyword argument. When set to
``True`` and ``applying_cuts`` is ``False``, the column callable is invoked only on
the events that pass all currently queued cuts (i.e. on ``selector.final``). The result is
then re-expanded back to the full data array size using masking, so the column remains
aligned with the unfiltered event array:

.. code-block:: python

    selector.set_column("MyDerivedQuantity", my_function, all_cuts=True)

This is useful when the column computation is expensive or only meaningful for events
that would survive the full selection, while still keeping the full array intact.
For example, in a dileptonic :math:`t\bar{t}` analysis, the reconstruction method only has 
:math:`90\%` efficiency but it may be of interest to also sometimes consider events where
the reconstruction method fails.

""""""""""""""""""""""""""""""""""""""""
When to use deferred cuts
""""""""""""""""""""""""""""""""""""""""

Deferring cut application is useful when you want to record the effect of multiple cuts
without committing to a particular order of removal - for example, when building a
complete cutflow that evaluates all cut combinations simultaneously, or when a derived
column must be computed before any rows are dropped but only needs to run on passing
events. The accumulated cuts and their combined boolean product remain available via
``selector.unapplied_cuts`` and ``selector.unapplied_product`` respectively.


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Ready to Run?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Great stuff! Check out the detailed guide on how on this page: :ref:`running-your-processor`.
