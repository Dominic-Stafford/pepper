.. _extending-config-example:

Extending the ``Config`` Class
==============================

JSON configuration is enough for most analyses, but some inputs do not fit naturally 
into a JSON file: ROOT files that need to be opened, objects that require non-trivial Python construction, or per-key validation that goes
beyond "is this a string". For these cases, Pepper's :class:`pepper.Config` class is designed to be subclassed.

This example walks through what inheritance/subclassing actually looks like by dissecting :class:`pepper.ConfigTTbarLL`, 
the configuration class shipped with Pepper for the dileptonic :math:`t\bar{t}` example analysis. 
For the broader picture of how the ``Config`` object is structured, see :doc:`/user_guide/config_object`.

.. note::
    Users are recomended to extend :class:`pepper.ConfigBasicPhysics` rather than implementing
    changes directly in the class. This ensures that merging with master branch updates is as 
    painless as possible, and that you get all the standard physics behaviours for free.


When to subclass
----------------

There are three main things that motivate writing your own ``Config`` subclass:

- **A new key whose value is not directly JSON-representable.** For example, the value the rest of your code wants 
  is a ROOT histogram or a ``correctionlib`` evaluator, but in the configuration file you would like to write just a path and a histogram name.
- **Per-key transformation.** You want ``config["my_key"]`` to return a parsed, ready-to-use object rather than the raw string or list that was
  written in the JSON file - and you want this transformation to happen exactly once, with the result cached.
- **A new required key.** You want Pepper to fail loudly at start-up if a required field is missing, rather than partway through processing.

If all you need is a new optional scalar or list of strings, you do not need to subclass - just read it directly from the config dictionary. 

.. note::

   Anything you can do with a custom ``Config`` subclass, you can also do with a plain JSON file and some helper 
   functions in your processor. The subclass just provides a convenient, centralised place to put that logic, 
   and a clean way to expose the processed values to the rest of your code.


Where each piece of behaviour lives
-----------------------------------

The ``Config`` class exposes three extension points, all of them simple dictionaries you mutate in ``__init__``:

``self.behaviors``
   Maps a config key name to a callable. When the key is accessed via ``config[key]``, the callable receives the 
   raw JSON value and returns the processed value that the rest of your code sees. Results are cached, so the
   callable runs at most once per key.

``self.special_vars``
   Maps a placeholder string (by convention prefixed with ``$``) to the name of another config key. When that 
   placeholder appears anywhere inside a string value, it is substituted at access time. The built-in ``$CONFDIR``,
   ``$STOREDIR``, and ``$DATADIR`` are defined here. Refer to :doc:`/user_guide/config_object` for a description
   of the built-in special variables.

``self.required_args``
   A list of key names that must be present. ``check_required_args`` will raise a :class:`pepper.config.ConfigError` if any are missing.

A subclass typically sets these up in ``__init__`` after calling ``super().__init__(...)``, so it inherits the parent's behaviours and adds to them.


A complete subclass example: ``ConfigTTbarLL``
-------------------------------------------------

We consider the configuration class ``ConfigTTbarLL`` shipped with Pepper. The full source is in 
:repo:`pepper/config_ttbarll.py`. The interesting part is the following snippet:

.. code-block:: python

   import numpy as np
   import uproot
   import hjson

   import pepper
   from pepper.scale_factors import ScaleFactors


   class ConfigTTbarLL(pepper.ConfigBasicPhysics):
       def __init__(self, path_or_file, textparser=hjson.load, cwd="."):
           super().__init__(path_or_file, textparser, cwd)

           self.behaviors.update({
               "drellyan_sf": self._get_drellyan_sf,
               "trigger_sfs": self._get_trigger_sfs,
           })

       def _get_drellyan_sf(self, value):
           if isinstance(value, list):
               path, histname = value
               with uproot.open(self._get_path(path)) as f:
                   hist = f[histname]
               return ScaleFactors.from_hist(hist)
           data = self._get_maybe_external(value)
           return ScaleFactors(
               bins=data["bins"],
               factors=np.array(data["factors"]),
               factors_up=np.array(data["factors_up"]),
               factors_down=np.array(data["factors_down"]),
           )

       def _get_trigger_sfs(self, value):
           path, histnames = value
           if len(histnames) != 3:
               raise pepper.config.ConfigError(
                   f"Need 3 histograms for trigger scale factors. Got {len(histnames)}"
               )
           ret = {}
           with uproot.open(self._get_path(path)) as f:
               for chn, histname in zip(("is_ee", "is_em", "is_mm"), histnames):
                   ret[chn] = ScaleFactors.from_hist(
                       f[histname], dimlabels=["lep1_pt", "lep2_pt"]
                   )
           return ret


There are a few things to note here:

- The class extends :class:`pepper.ConfigBasicPhysics` rather than the bare :class:`pepper.Config`. 
  This is the right base class for any analysis using NanoAOD - it adds all the standard physics 
  behaviours (jet corrections, electron/muon SFs, b-tag SFs, etc.) on top of the base.
  Only extend :class:`pepper.Config` directly if you are deliberately bypassing those standard behaviours.
- ``super().__init__`` runs first, so by the time the subclass touches ``self.behaviors`` it already 
  contains every entry inherited from the parent. ``self.behaviors.update({...})`` adds the new keys on top without
  losing anything.
- The behaviour callables use the inherited helpers ``self._get_path`` and ``self._get_maybe_external``. 
  ``_get_path`` resolves a path relative to the config file's ``cwd``, expanding ``$CONFDIR`` and friends.
  ``_get_maybe_external`` accepts either an inline value or a path to a separate JSON file, and parses it transparently. 
  Reusing these helpers keeps your subclass consistent with how the rest of Pepper handles paths and external files.
- Errors are raised as :class:`pepper.config.ConfigError`, not ``ValueError`` or ``RuntimeError``. 
  Pepper formats those errors nicely at the top level, so users see the configuration problem rather than a Python traceback.

The corresponding configuration JSON snippet looks like:

.. code-block:: json

   {
       "drellyan_sf": ["$CONFDIR/dy_sf.root", "DY_SF_2018"],
       "trigger_sfs": [
           "$CONFDIR/trigger_sf.root",
           ["sf_ee", "sf_emu", "sf_mm"]
       ]
   }

When the processor accesses ``config["drellyan_sf"]``, it receives a fully constructed 
:class:`pepper.scale_factors.ScaleFactors` object - never the raw list - and the ROOT file is opened only once. 


.. note::
    Config behaviours are lazy: if the processor never accesses ``config["trigger_sfs"]``, 
    the ROOT file containing the trigger scale factors is never opened, and the histograms are never read. 
    This can save time and memory if you have a large number of expensive-to-construct config entries 
    but only use a subset of them in a given run.


.. note::
    Config lookups are cached: if the processor accesses ``config["drellyan_sf"]`` multiple times, 
    the ROOT file is opened only once and the same ``ScaleFactors`` object is returned each time.



Adding a required key and a special variable
---------------------------------------------

The ``ConfigTTbarLL`` example only adds behaviours. Suppose you also want to require a new key called 
``analysis_region`` and to introduce a new placeholder ``$REGIONDIR`` that expands to a directory 
derived from a new key, say ``region_inputs_dir``. 
You can do both of these in the same subclass by mutating ``self.required_args`` and ``self.special_vars``:

.. code-block:: python

   class ConfigMyAnalysis(pepper.ConfigBasicPhysics):
       def __init__(self, path_or_file, textparser=hjson.load, cwd="."):
           super().__init__(path_or_file, textparser, cwd)

           self.required_args.append("analysis_region")

           self.special_vars["$REGIONDIR"] = "region_inputs_dir"

           self.behaviors.update({
               "my_input_table": self._get_maybe_external,
           })

           self.check_required_args()

The call to ``check_required_args`` at the end of ``__init__`` ensures that missing keys fail fast at construction time, before any analysis code runs.
The parent class does not call this for you - it is your responsibility to trigger the check after you have finished extending ``required_args``.


Pointing the processor at your subclass
---------------------------------------

A processor selects its config class by setting ``config_class`` at the class
level:

.. code-block:: python

   class MyProcessor(pepper.ProcessorBasicPhysics):
       config_class = ConfigMyAnalysis

       def process_selection(self, selector, dsname, is_mc, filler):
           ...

When ``pepper.runproc`` loads your processor and configuration JSON, it will instantiate ``ConfigMyAnalysis`` 
rather than the inherited ``ConfigBasicPhysics``, and the new keys, behaviours, and special variables
become available to your selection code.


See also
--------

- :doc:`/user_guide/config_object` - conceptual description of the ``Config`` class.
- :ref:`configuration-reference` - the configuration keys recognised by ``Config`` and ``ConfigBasicPhysics`` out of the box.
- :ref:`config-inheritance-example` - reuse JSON configurations across eras or analyses without writing any Python.
- :repo:`pepper/config_ttbarll.py` - the source for the example dissected above.
- :repo:`pepper/config_basic.py` - a much larger subclass that defines all the standard HEP behaviours, 
  useful as a reference when writing involved custom behaviours of your own.