.. _bundled-scripts-example:

Using the Bundled Scripts
=========================

A complete Pepper analysis is rarely just a single call to ``pepper.runproc`` with a processor.
Several pre- and post-processing steps are needed to produce auxiliary inputs (lumi factors, pileup weights, b-tagging efficiencies),
to debug the implementation of corrections in the main processor, to produce plots from the output histograms, 
and to inspect the outputs (plotting, file management, cluster housekeeping). All of these are shipped as scripts under :repo:`scripts/`.

The reference :doc:`/user_guide/scripts` lists every script bundled in Pepper with a one-line description. 
This page instead walks through them in the order you would actually use them when starting a new analysis 
from scratch, explaining how they fit together.


The big picture
---------------

A typical analysis loop looks roughly like this:

.. code-block:: text

   1.  Set up datasets + cross sections
   2.  compute_mc_lumifactors.py        -> lumifactors.json   (or use "mc_lumifactors": false)
   3.  compute_pileup_weights.py        -> pileup.root        (Analyses without correction-
                                                               lib pile up correction only)
   4.  Run main processor (debug)       -> sanity check
   5.  Run main processor (full, no btag_sf)
                                        -> btageff histogram  (if b-tagging)
   6.  generate_btag_efficiencies.py    -> btag_eff.root      (if b-tagging)
   7.  Run main processor (full)        -> full histograms
   8.  make_plots.py                    -> publication plots
   9.  delete_duplicate_outputs.py      -> cleanup if needed

Steps 2, 3, and 6 produce inputs that are then consumed by the main processor in step 7. 
They typically only need to be re-run when the input datasets or POG recommendations change.


Step 1: Datasets and cross sections
-----------------------------------

Two JSON files form the basis of every Pepper analysis:

- ``crosssections.json`` - a flat ``{dataset_name: cross_section_pb}`` dictionary. See :repo:`example/crosssections.json`.
- The main config file's ``exp_datasets`` and ``mc_datasets`` blocks - the list of CMS datasets to run over. See :repo:`example/example_config.json`.

If your data lives on a remote site, also configure XRootD access at this point (see :doc:`/user_guide/remote_data`). 
For analyses that will run on a lot of data, requesting a local Rucio transfer up front can save you a great deal of time downstream:

.. code-block:: bash

   python scripts/rucio_create_rules.py my_config.json T2_DE_DESY

Rucio rules need to be manually approved by the site after creation, so this step should be done well in advance of when you need the data.
Expect the transfer to take a few working days to be approved and completed.


Step 2: Compute MC lumi factors
-------------------------------

For a detailed discussion of MC normalisation in Pepper, see :doc:`/user_guide/mc_normalization`. There are two main 
approaches to MC normalisation in Pepper:

1. The ``"mc_lumifactors": false`` approach, where the lumifactors are computed on the fly by the processor. **This 
   is the recommended approach for first time users**. If using this, you can skip step 2 entirely and just set 
   ``"mc_lumifactors": false`` in your config.
#. The "lumifactors" approach, where you compute the :math:`\mathcal{L}\sigma / \sum w_{\mathrm{gen}}` factor for each 
   MC sample and apply it as a single weight in the processor. This requires that you compute the factors yourself using 
   ``compute_mc_lumifactors.py``.

``compute_mc_lumifactors.py`` produces the :math:`\mathcal{L}\sigma / \sum w_{\mathrm{gen}}` factor for each MC
sample, which is what scales the simulated yields to match the data luminosity.

.. code-block:: bash

   python scripts/compute_mc_lumifactors.py my_config.json lumifactors.json

The script reads the MC dataset list and luminosity from your config, opens each NanoAOD file to read the sum of 
generator weights, and writes the combined factor for each dataset. You can submit this on HTCondor with ``--condor N`` 
to run multiple datasets in parallel.

The output goes into your config as:

.. code-block:: json

   "mc_lumifactors": "$CONFDIR/lumifactors.json"

Re-run this whenever you add or remove MC datasets, change cross sections, or change the integrated luminosity.

.. tip::

   Pass ``--skip`` if you have already computed factors for some datasets and only want to fill in 
   new ones. This is much faster than starting from scratch.


Step 3: Pileup weights (Run 2)
------------------------------

For Run 2 analyses that do not use the correctionlib-based pileup correction, ``compute_pileup_weights.py`` computes data/MC pileup ratios
from MC pileup distributions and externally provided data pileup histograms.

The script is itself a Pepper processor and is invoked with the same machinery as your main processor:

.. code-block:: bash

   python -m pepper.runproc scripts/compute_pileup_weights.py my_config.json \
       --output pileup_output

Your config must declare ``data_pu_hist``, ``data_pu_hist_up``, and ``data_pu_hist_down`` pointing to the data pileup ROOT files. The output
directory will contain ``pileup.root``, which then plugs into your config:

.. code-block:: json

   "pileup_reweighting": "$CONFDIR/pileup.root"

Newer analyses can skip this script entirely and use the correctionlib (if available) approach described in :doc:`/user_guide/scale_factors`. 


Step 4: First test run of the main processor
---------------------------------------------

Before producing any analysis-specific inputs, run the main processor in debug mode with the configuration assembled so far:

.. code-block:: bash

   python -m pepper.runproc my_processor.py my_config.json --debug --output test_out

``--debug`` processes only one chunk per dataset, which is enough to verify that all paths resolve and all corrections load and are applied.


Step 5: First production pass (for b-tagging only)
--------------------------------------------------
 
Skip this step if your analysis does not use b-tagging.
 
B-tagging scale factors require per-analysis (phase-space dependent) efficiency histograms measured with 
full statistics on your actual selection. To produce them, run the main processor over the full sample with 
the ``btageff`` histogram enabled in ``hist_config.json`` (see :repo:`example/hist_config.json` for the definition) 
and **without** ``btag_sf`` in the config:
 
.. code-block:: bash
 
   python -m pepper.runproc my_processor.py my_config.json \
       --condor 200 --output btageff_run
 
The analysis histograms produced by this run will be missing b-tagging weights and should be discarded - only the ``btageff`` histogram is needed.
for the next step.
 
 
Step 6: B-tagging efficiencies
------------------------------
 
Skip this step if your analysis does not use b-tagging.
 
Hand the ``hists.json`` from step 5 to ``generate_btag_efficiencies.py``:
 
.. code-block:: bash
 
   python scripts/generate_btag_efficiencies.py btageff_run/hists.json btag_eff.root \
       --cut HasJets
 
The ``--cut`` argument names the cut applied immediately *before* the b-tagging requirement in your selection.
The default cut before the b-tagging requirement is ``HasJets``. Thescript reads the ``btageff`` histogram and writes a small ROOT file with
central efficiencies and systematic variations. You can skip the variations if you only want to produce central b-tagging SFs by passing ``--central``.
 
Plug the result into your config:
 
.. code-block:: json
 
   "btag_sf": [["/cvmfs/.../btagging.json.gz", "$CONFDIR/btag_eff.root"]]
 
For the full end-to-end SF setup, see :ref:`scale-factors-example`.
 
 
Step 7: Production run
----------------------
 
With all inputs in place, submit the full processing job. For anything
beyond a few datasets, this means HTCondor:
 
.. code-block:: bash
 
   python -m pepper.runproc my_processor.py my_config.json \
       --condor 200 --output prod_v1
 
A few practical suggestions for production runs:
 
- Use ``tmux`` or ``nohup`` so the controlling process survives a closed terminal. 
  ``pepper.runproc`` itself stays alive until all Condor jobsfinish.
- Keep an eye on ``pepper_logs/``.
- If a run fails partway through, **you do not need to start over**. Re-run the same command with ``--resume`` 
  and Pepper will pick up where it left off (see :ref:`running-your-processor`). However, this requires that
  you also run ``delete_duplicate_outputs.py`` after the production finishes to remove any output files that 
  were produced more than once. See Step 9 for details.
 
 
Step 8: Plotting
----------------
 
``make_plots.py`` is the standard way to visualise the output:
 
.. code-block:: bash
 
   python scripts/make_plots.py example_plotting.json prod_v1/hists.json --outdir plots
 
The plotting configuration (:repo:`example/example_plotting.json` is a small starting point) groups MC datasets 
into stacked categories, defines data, controls rebinning, axis labels, and signal overlays. The script reads each 
histogram from ``hists.json`` and writes one plot per histogram into the output directory.
 
Pass ``--log`` for log-scale axes and ``--ext png`` (or ``svg``, ``pdf``) to control the output format. 
``--cut`` restricts plotting to a single cutflow stage, which is useful for iterating on plot styling without
producing the full set every time.
 
The script is deliberately not a one-size-fits-all solution. Once you need plots that fall outside its capabilities, 
the recommended path is to copy the script into your analysis repository and adapt it - the helper functions in 
:mod:`pepper.hist_utils` (rebinning, grouping, projecting) are intended for direct use in your own plotting code as well.
 
 
Step 9: Cleanup
---------------
 
Two scripts handle common housekeeping needs:
 
- ``delete_duplicate_outputs.py`` - identify and remove per-event output files that were produced more than once. 
  This commonly happens after a ``--resume`` run picks up jobs that had partially completed before the failure. 
  Run this once after the production batch finishes to avoid double-counting.
- ``get_bad_local_files.py`` - find NanoAOD files in your local store that exist on disk but cannot be opened 
  (corrupted transfers, permission problems, etc.). Useful to run periodically on a long-lived store directory.
 
 
Analysis-specific scripts
-------------------------
 
The scripts directory also contains several files prefixed with ``ttbarll_`` that are specific to the 
:math:`t\bar{t}\to\ell\ell` reference analysis shipped with Pepper -- DY scale factor production, kinematic
reconstruction histograms. Treat these as worked templates: if your analysis needs an equivalent, 
copy the relevant script into your own analysis repository and adapt it. They are not designed to be invoked unmodified from a different analysis.

Other Scripts
^^^^^^^^^^^^^
Some scripts are not stored in the main repository but are instead found in separate repositories for specific analyses. 
The following list, which is not necessarily exhaustive, 
includes scripts that are not general-purpose enough to be included in the main repo but are still useful as examples for how to implement common analysis tasks:

- Cross trigger scale factor production: :repo-gitlab:`https://gitlab.cern.ch/cms-analysis/general/pepper/examples/cross-trigger-sfs`
- Assemble a pdf of plots outputted by Pepper: :repo-gitlab:`https://gitlab.cern.ch/mabaattr/pepper-output-plotting`


See also
--------
 
- :doc:`/user_guide/scripts` - one-line reference for every script.
- :doc:`/user_guide/running_a_processor` - in-depth coverage of
  ``pepper.runproc`` and HTCondor submission.
- :ref:`scale-factors-example` - the workflow above, focused on getting
  the corrections themselves right.
