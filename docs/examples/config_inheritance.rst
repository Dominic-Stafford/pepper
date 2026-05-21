.. _config-inheritance-example:

Config Inheritance
==================

Pepper configuration files support inheritance through the ``import`` key.
This lets you factor settings that are common to several analyses, or to
several data-taking eras of the same analysis, into a single base file and
keep the per-analysis or per-era files as small overrides.

This example walks through a typical multi-era setup. For a more general
discussion of the ``Config`` object and how to extend it from Python, see
:doc:`/user_guide/config_object`.


When to use ``import``
----------------------

Reach for config inheritance whenever the same JSON keys would otherwise be
duplicated across several configuration files. Common cases include:

- A base ``ttbar`` selection that is shared between several signal regions or
  control regions.
- A multi-era analysis where most of the selection is era-independent, but
  luminosity, golden JSON, trigger lists, scale-factor inputs, and dataset
  lists differ between 2016, 2017, and 2018.
- A team-wide "house defaults" file (object kinematic cuts, MET filters,
  plotting style) that individual analyses extend.

If the only thing varying between two configurations is which datasets are
processed, you usually do not need ``import`` - you can just point
``pepper.runproc`` at a different ``--datasets`` file. ``import`` is most
useful when many keys differ together.


A two-file example
------------------

The base file collects everything that does not depend on the data-taking
year. The era-specific file imports it and overrides only the keys that
actually change.

``base.json`` -- shared settings:

.. code-block:: json

   {
       "datadir": ".",
       "store": "/pnfs/desy.de/cms/tier2/store",

       "ele_eta_min": -2.4,
       "ele_eta_max": 2.4,
       "good_ele_pt_min": 20.0,

       "muon_eta_min": -2.4,
       "muon_eta_max": 2.4,
       "good_muon_pt_min": 20.0,

       "compute_systematics": false,
       "hists": "$CONFDIR/example_histograms.json"
   }

``config_2018.json`` -- era-specific overrides:

.. code-block:: json

   {
       "import": "/abs/path/to/base.json",

       "year": "2018",
       "luminosity": 59.74,
       "lumimask": "$CONFDIR/Cert_314472-325175_..._JSON.txt",

       "exp_datasets": {
           "MuonEG":      ["/MuonEG/Run2018A-02Apr2020-v1/NANOAOD",  "..."],
           "DoubleMuon":  ["/DoubleMuon/Run2018A-02Apr2020-v1/NANOAOD", "..."],
           "EGamma":      ["/EGamma/Run2018A-02Apr2020-v1/NANOAOD",  "..."]
       },
       "mc_datasets": {
           "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8": ["/TTTo2L2Nu_..._v1/NANOAODSIM"]
       },
       "mc_lumifactors": false
   }

Pointing ``pepper.runproc`` at ``config_2018.json`` will load ``base.json``
first, then apply everything in ``config_2018.json`` on top of it. A second
file ``config_2017.json`` would have the same shape with different values for
``year``, ``luminosity``, ``lumimask``, and the dataset lists.

.. hint::

   The ``import`` key is processed *before* the special-variable substitution
   that turns ``$CONFDIR`` into the directory of the configuration file.
   ``$CONFDIR`` therefore cannot be used inside the value of ``import``
   itself - use either an absolute path, or a path relative to the directory
   from which you launch ``pepper.runproc``.


How values are merged
---------------------

The merge is a single-level (shallow) update. Keys present in the importing
file replace the corresponding keys in the imported file outright; keys not
present in the importing file are inherited unchanged.

This means nested objects are **not** deep-merged. If you override a key whose
value is itself a dictionary, you must repeat any nested entries you want to
keep.

For example, suppose ``base.json`` defines:

.. code-block:: json

   {
       "dataset_trigger_map": {
           "MuonEG":     ["HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL"],
           "Muon":       ["HLT_IsoMu30"],
           "EGamma":     ["HLT_Ele27_WPTight_Gsf"]
       }
   }

and the era file overrides only the ``EGamma`` entry:

.. code-block:: json

   {
       "import": "/abs/path/to/base.json",
       "dataset_trigger_map": {
           "EGamma": ["HLT_Ele28_WPTight_Gsf"]
       }
   }

The merged ``dataset_trigger_map`` will contain only the ``EGamma`` key.
The ``MuonEG`` and ``Muon`` entries from ``base.json`` are lost,
because the whole ``dataset_trigger_map`` value was replaced. To keep them,
copy them into the era file as well.

.. note::

   This behaviour is intentional: it keeps the merge semantics simple and
   predictable. If you find yourself copying large nested objects back and
   forth, consider moving them into a separate JSON file and referencing
   that file by path instead.


Chained and circular imports
----------------------------

An imported file may itself contain an ``import`` key, so you can build chains
of arbitrary length:

.. code-block::

   config_2018_signal.json -> base_2018.json -> base_common.json

At each level, the same shallow-update rule applies: the file closer to the
one you actually load wins.

Circular imports are detected and rejected with a ``ConfigError``, so you
cannot accidentally create an infinite loop by having two files import each
other.


See also
--------

- :ref:`configuration-reference` - the full list of configuration keys.
- :doc:`/user_guide/config_object` - the ``Config`` class itself, including
  how to subclass it to add fields that ``import`` cannot express.
- :repo:`example/example_config.json` - the single-file example shipped with
  Pepper, useful as a starting point for the base file in your own setup.