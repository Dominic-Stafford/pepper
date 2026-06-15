Memory Optimization
========================================
In general, one should try to run with a large chunksize (i.e. number of events processed per chunk) for optimal performance.
The default in pepper is 500000. However, resident memory consumption of course increases with the chunksize. In case you run
into issues with condor workers being killed for running out of memory, you can try several solutions:

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Histogram optimization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Large histograms can take up a lot of memory especially when running with many systematic uncertainties and/or many categories.
To reduce the number of histograms, you can

- fill histograms only for certain cuts (e.g. at the end of the selection) using ``cuts_to_histogram`` in the config,
- consider systematic uncertainties only for certain histograms (e.g. the ones used for a fit) by specifying ``do_systs: false`` in the histogram definition,
- fill certain histograms only for certain categories by specifying ``only_for_cats: ...``.

For all these options, see :ref:`configuration-reference`.

.. _column-unloading:

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Freeing Memory with Column Unloading
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, pepper will keep all data that it loads from nanoAOD in memory
for the entire processing of the chunk. This can be wasteful if one does not need the raw data in
the analysis, but only objects derived from it. For example, when working with the very large
``GenPart`` column, one frequently only needs a few specific particles, e.g. gen-level top quarks,
which one then will set as a separate column (see method ``gentop`` in :repo:`pepper/processor_basic.py`.)
The original ``GenPart`` column is then no longer needed, but will by default still be kept in memory.

To fix this, the ``Processor.unload_column`` method removes a column and its count branch from the cache:

.. code-block:: python

    # After the GenPart-based top-quark column has been built, free the raw branches
    self.unload_column("GenPart")

    # After the combined Lepton column is available, free the raw object branches
    self.unload_column("Electron")
    self.unload_column("Muon")

Matching is greedy: ``unload_column("Jet")`` removes every cached key whose path contains
``"/Jet"`` or ``"/nJet"``, so it covers all Jet sub-fields (``Jet_pt``, ``Jet_eta``, …) and
the count array ``nJet`` at once.

A good pattern is to call ``unload_column`` immediately after the last ``set_column`` or
``add_cut`` that depends on a given raw branch. The ``processor_ttbarll.py`` processor shows
this in practice:

.. code-block:: python

    # GenPart only needed for the generator-level top column:
    selector.set_column("gent_lc", self.gentop)
    self.unload_column("GenPart")

    # Electron/Muon only needed until the combined Lepton column is set:
    selector.add_cut("AtLeast2Leps", partial(self.lepton_pair, is_mc), no_callback=True)
    self.unload_column("Electron")
    self.unload_column("Muon")

    # Jet/MET columns are reused for each JEC variation, so unload only after
    # the nominal variation is complete:
    if variation.name is None:
        self.unload_column("Jet")
        self.unload_column("GenJet")
        self.unload_column("MET")

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Reducing the chunksize
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If none of these help, you can reduce the chunksize using the ``--chunksize`` command line argument. Note, however,
that this usually comes at a significant decrease in performance due to additional overhead.