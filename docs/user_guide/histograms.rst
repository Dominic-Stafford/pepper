.. _histograms:

Histograms
==========

Histograms are Pepper's primary form of analysis output. They are defined declaratively in JSON. 
They are written once per processor run, and read back with a small set of helper utilities. 
This page covers the full lifecycle: how histograms are defined, what files Pepper produces, 
how to load and manipulate them, and how to keep the output size manageable.

The exhaustive list of configuration keys is in the :ref:`configuration-reference` 
-- in particular the *Histogram Definition* section. This page assumes you have at least 
skimmed the introductory paragraphs in :doc:`conceptual_overview` and focuses on the *how*.

Defining histograms
-------------------

Histograms are listed in the JSON file referenced by the ``hists`` config key and it is 
typically a separate file like ``hist_config.json``. 
Each entry has a name (the dictionary key), one or more axis definitions in ``bins``, 
and a ``fill`` block that says where to get the values:

.. code-block:: json

   "leading_lepton_pt": {
       "bins": [
           {
               "name": "pt",
               "label": "Lepton $p_{\\mathrm{T}}$",
               "n_or_arr": 100,
               "lo": 0,
               "hi": 400,
               "unit": "GeV"
           }
       ],
       "fill": {
           "pt": ["Lepton", "pt", {"leading": 1}]
       }
   }

Pepper will produce one histogram per cut for which the ``fill`` data is defined, 
with implicit category axes added automatically.

The fill list is a *DataPicker*. It walks an awkward-array starting from the per-event record: 
each string item becomes attribute or key access, and each dictionary item invokes a named operation. 
The example above reads ``events.Lepton.pt``, then keeps only the leading entry per event.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Common DataPicker patterns
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A handful of patterns cover most needs:

**Per-object kinematics.**
   ``["Jet", "pt", {"leading": 2}]`` -- the second-leading jet's
   :math:`p_T`. ``{"leading": [1, N]}`` keeps the highest-:math:`p_T` ``N``
   objects; combine with a sum or other aggregation if you want one value
   per event.

**Counting.**
   ``["Jet", {"function": "num"}]`` returns the per-event jet
   multiplicity. ``num`` and ``sum`` both work on ragged arrays.

**Element-wise functions.**
   ``["MET", "pt", {"function": "log10"}]``. The full list of available
   functions is in :mod:`pepper.hist_defns` -- ``sin``, ``cos``, ``abs``,
   ``log``, ``sqrt``, ``int``, etc., plus a few HEP-specific helpers like
   ``leaddiff`` (difference between the two leading entries).

**Multi-dimensional histograms.**
   Add more entries to ``bins`` and one ``fill`` entry per axis:

   .. code-block:: json

      "jet_pt_vs_eta": {
          "bins": [
              {"name": "pt",  "label": "Jet $p_T$ (GeV)", "n_or_arr": 50, "lo": 0, "hi": 500},
              {"name": "eta", "label": "Jet $\\eta$",     "n_or_arr": 50, "lo": -2.5, "hi": 2.5}
          ],
          "fill": {
              "pt":  ["Jet", "pt"],
              "eta": ["Jet", "eta"]
          }
      }

**Categorical splits.**
   ``cats`` defines extra string-category axes. Each entry in ``fill`` for a category lists the named values 
   and the boolean DataPicker that selects events in that value. The b-tagging efficiency histogram in
   :repo:`example/hist_config.json` is the canonical example: a single ``btagged`` axis with ``"yes"`` 
   and ``"no"`` keys, filled from the jet's ``btagged`` field and its negation.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Restricting when histograms fill
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, every histogram is produced after every cut whose state is sufficient to evaluate its fill. 
This is convenient during analysis development but produces a large number of files for mature analyses. Two keys narrow this:

``step_requirement`` (per-histogram)
   Restrict a histogram to fill only after a named cut or once a named column exists. Use ``"cut:HasJets"`` to 
   defer until ``HasJets`` has run, or ``"column:reco_top"`` to defer until that column has been computed.

``cuts_to_histogram`` (top-level)
   A list of cut names. If set, histograms are produced *only* at those cuts, regardless of what individual histograms 
   would otherwise allow. This is the more aggressive lever and the one to reach for when output size becomes a problem.

Custom weights are supported via ``weight`` (replaces the event weight entirely) and ``weight_factor`` 
(multiplies the event weight, including its systematic variations). Use ``weight_factor`` for analysis-specific
reweightings that should preserve the standard systematics on top.

Output layout
-------------

After ``pepper.runproc`` finishes, the output directory contains:

.. code-block:: text

   output/
   |-- cutflows.json              # weighted/unweighted yields per cut
   |-- hists/
   |   |-- Cut_000_NoCut_leading_lepton_pt.coffea
   |   |-- Cut_001_HasLeptons_leading_lepton_pt.coffea
   |   |-- Cut_001_HasLeptons_njet.coffea
   |   |-- ...
   |-- ...

Histogram files are named ``Cut_NNN_<cutname>_<histname>.<ext>``. The ``NNN`` prefix preserves cut order numerically, 
since dictionary key order is not guaranteed across implementations.

The default extension is ``.coffea``. Two configuration choices change this:

``"histogram_format": "hist"``  (default)
   Each histogram is a serialised :mod:`hist` object inside a coffea container. 
   Fastest to write, fastest to load back into Python, no external tooling required. 
   This is the default and the recommended format for most users.

``"histogram_format": "coffea"``
   Same on-disk format as ``hist`` but kept as a label for forward compatibility. Identical to ``hist``.

``"histogram_format": "root"``
   Writes ``.root`` files instead. Each histogram is split across multiple ``TH1`` / ``TH2`` objects, 
   one per combination of category-axis values, since ROOT does not support string-category axes natively. 
   Use this when downstream tooling outside Pepper expects ROOT input - otherwise just keep the default.

Reading histograms back
-----------------------

The :class:`pepper.HistCollection` class is a thin ``Mapping`` over the ``hists.json`` index. 
It does not load any histogram data until you explicitly ask for it -- iterating over the keys is cheap, 
even for runs with many histograms.

.. code-block:: python

   import pepper

   with open("output/hists/hists.json") as f:
       hc = pepper.HistCollection.from_json(f)

   # Keys are tuples of (cut_name, histogram_name)
   for cut, hist_name in hc.keys():
       print(cut, hist_name)

   # Load one histogram on demand
   h = hc.load(("HasLeptons", "leading_lepton_pt"))

Loading a single file directly, without going through the index, is also
supported:

.. code-block:: python

   h = pepper.HistCollection.from_single_hist(
       "output/hists/Cut_001_HasLeptons_leading_lepton_pt.coffea"
   )

The loaded ``h`` is a :class:`hist.Hist` object. Project, slice, and plot it using the 
standard :mod:`hist` and :mod:`mplhep` interfaces -- Pepper does not wrap them.

For ROOT-format output, ``HistCollection`` reassembles the per-category sub-histograms into a 
single ``hist.Hist`` on load, so user code does not need different paths for the two formats.

Manipulating histograms
-----------------------

The :mod:`pepper.hist_utils` module collects four operations that the :mod:`hist` library either does 
not provide or provides in a less convenient form. All four return a new histogram 
(except ``replace_missing_systematics`` and ``scale_histogram``, which operate in place).
None of them is specific to Pepper output and they can be used on any ``hist.Hist`` with the expected axes.

``group_histogram(h, old_axis, new_axis, mapping)``
   Merge categories on ``old_axis`` into broader groups on a newly created ``new_axis``. 
   The typical use is collapsing the ``dataset`` axis (one MC sample per category) into a process axis 
   (e.g. ``"ttbar"``, ``"single_top"``, ``"DY"``, ``"diboson"``). 
   The mapping is ``{new_category: [list_of_old_categories]}``.

``rebin_histogram(h, axis_name, edges, new_axis_label=None)``
   Rebin a single dense axis to a coarser binning. ``edges`` may be an integer (an :mod:`hist` ``rebin`` factor) 
   or an explicit array of edges. Explicit edges must be a subset of the original edges 
   -- non-aligned rebinnings raise ``ValueError`` rather than silently interpolating.

``replace_missing_systematics(h)``
   Fixes a specific quirk of Pepper output: alternate-dataset systematics (top mass, ``hdamp``, tune, etc.) 
   are filled only on their target dataset and zero everywhere else, which produces wrong results when the 
   ``dataset`` axis is collapsed across processes. Call this *before* ``group_histogram`` whenever both 
   are needed. Modifies the histogram in place.

``scale_histogram(h, axis_name, scales)``
   Multiply specific category values on ``axis_name`` by per-category scale factors. ``scales`` is a 
   ``{category_value: factor}`` dictionary. Mostly useful for ad-hoc post-fit corrections during plotting iteration.

A typical post-processing chain looks like:

.. code-block:: python

   from pepper import HistCollection
   from pepper.hist_utils import (
       replace_missing_systematics, group_histogram, rebin_histogram,
   )

   with open("output/hists/hists.json") as f:
       hc = HistCollection.from_json(f)
   h = hc.load(("HasLeptons", "leading_lepton_pt"))

   replace_missing_systematics(h)
   h = group_histogram(h, "dataset", "process", {
       "ttbar":      ["TTTo2L2Nu", "TTToSemiLeptonic", "TTToHadronic"],
       "single_top": ["ST_tW_top", "ST_tW_antitop", "ST_t-channel_top"],
       "DY":         ["DYJetsToLL_M-10to50", "DYJetsToLL_M-50"],
   })
   h = rebin_histogram(h, "pt", [0, 30, 60, 100, 150, 250, 400])

The bundled ``make_plots.py`` script (see :ref:`bundled-scripts-example`) wraps a chain very similar to this 
one with JSON-driven configuration. For one-off plots and exploratory work, calling these helpers directly
in a notebook is usually faster than maintaining a plotting config.

See also
--------

- :ref:`configuration-reference` -- exhaustive list of histogram configuration keys.
- :doc:`conceptual_overview` -- where histograms fit in the overall Pepper architecture.
- :ref:`bundled-scripts-example` -- where ``hists.json`` fits into the full analysis workflow.
- :repo:`example/hist_config.json` -- a worked configuration covering most of the patterns described above, including categorical splits.