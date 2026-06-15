.. _scale-factors-example:

Applying Scale Factors
======================

This example walks through the practical workflow of adding scale factors to an analysis configuration. 
The reference page :doc:`/user_guide/scale_factors` documents the exact key syntax for every
scale factor Pepper supports -- this page focuses on the order of operations, 
where the input files come from, and how to verify the result.

We build up the configuration incrementally so that you can run the processor after each step 
and confirm that the new corrections behave as expected before adding the next one. 
The configuration snippets below are taken from Pepper's 2024 integration-test configuration with minor changes. 
The corresponding ``crosssections.json``, dataset blocks, selection cuts, and channel/
trigger maps are not shown here -- this page concentrates on the correction-related keys.

.. tip::

   Pepper configuration files are parsed with :mod:`hjson`, which
   accepts ``//`` and ``#`` comments inside JSON. The snippets below use
   this freely. Plain ``json.load`` will reject these files; use
   :func:`hjson.load` if you need to read them from your own scripts.

We will be referencing the POG-recommended correction files for 2024 Run 3 analyses, which are all available on CVMFS.
The easiest way to locate which corrections to use is the `CMS Analysis Corrections`_ website. 
It provides a user-friendly interface to browse the available corrections, find the right ones for your era and analysis type. 
Moreover, all corrections are version-tagged to ensure corrections are not updated without notice.

.. _CMS Analysis Corrections: https://cms-analysis-corrections.docs.cern.ch/

Before you start
----------------

You will need a working configuration that already runs without scale factors 
-- the :repo:`example/example_config.json` shipped with Pepper is a good starting point. 
Confirm that

.. code-block:: bash

   python -m pepper.runproc example_processor.py example_config.json --debug

completes without errors. ``--debug`` restricts the run to one chunk per dataset.

Throughout this example we assume a 2024 Run 3 analysis. The procedure is the same for other eras. 
We use the POG-recommended correctionlib files for each era as listed on the
`CMS Analysis Corrections`_ website. 
For Run 2 UL analyses, the procedure is similar but you will be using a mix of correctionlib and legacy ROOT files 
since not all Run 2 UL corrections are available in correctionlib format.

In the following, we assume a repository structure of the form:

.. code-block:: text
    
   event-selection/
   ├── config/  # This is $CONFDIR
   │   ├── example_config.json
   │   └── ...
   ├── data/  # This is $DATADIR
   │   ├── btag/
   │   │   └── btag_sf.root
   │   └── ...
   ├── pepper/
   │   ├── scripts/
   │   |   └── generate_btag_efficiencies.py
   │   └── ...
   └── ...

Step 1: Data eras, cross sections and luminosity
------------------------------------------------

Initially, the configuration needs to specify the data eras, the MC cross sections, and the integrated luminosity.



.. code-block:: json

   {
       "year": "2024",
       "luminosity": 109.95,
       "crosssections": "$CONFDIR/crosssections.json",
       "lumimask": "$CONFDIR/Cert_Collisions2024_378981_386951_Golden.json",
       "mc_lumifactors": false,
       "data_eras": {
           "2024B":   [378971, 379411],
           "2024C":   [379412, 380252],
           "2024D":   [380253, 380947],
           "2024Ev1": [380948, 381383],
           "2024Ev2": [381384, 381943],
           "2024F":   [381944, 383779],
           "2024G":   [383780, 385813],
           "2024H":   [385814, 386408],
           "2024Iv1": [386409, 386797],
           "2024Iv2": [386798, 387121]
       }
   }

``crosssections.json`` is a flat dictionary mapping MC dataset names to cross sections in fb. 
:repo:`example/crosssections.json` shows the format. We set ``mc_lumifactors`` to ``false`` 
to compute the lumi-factors at runtime.

The ``lumimask`` key points to the JSON file defining the certified luminosity for the data eras. This can
for example be found on the `PdmV Run 3 TWiki <https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun3Analysis>`_. 
It defines which luminosity sections are good for analysis and which are not.

The ``data_eras`` block maps each data era to a list of run ranges.
This is needed for the processor to determine which era-specific corrections to apply to each event.


Step 2: Pileup reweighting
--------------------------

Pileup reweighting is the next universal correction. The Run 3 path is correctionlib-based:
We find the latest pileup weights on the `CMS Analysis Corrections`_ website, under the LUM POG section. 
The configuration looks like:

.. code-block:: json

   {
       "pileup_reweighting": [
           "/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/2026-04-15/puWeights_BCDEFGHI.json.gz",
           "Collisions24_BCDEFGHI_goldenJSON"
       ]
   }

The exact correction name depends on the era. After adding this and re-running with ``--debug``, the cutflow
should now report a non-trivial weighted yield even for events that pass no cuts, because every MC event carries a pileup weight.

If you are working on a Run 2 analysis without correctionlib pileup files, \
use the legacy ROOT-file approach documented in the *Pile-Up* section of :doc:`/user_guide/scale_factors`.


Step 3: Lepton ID, isolation, and reconstruction
------------------------------------------------

Lepton corrections are configured as lists of corrections that are multiplied together. 
Each entry can refer to a correctionlib JSON or (for older eras) a ROOT file. A typical 2024 dilepton setup looks like:

.. code-block:: json

   {
        "electron_sf": [
            [
                "/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/2025-12-15/electron.json.gz",
                "Electron-ID-SF",
                {"year": "2024Prompt", "WorkingPoint": "wp90iso"}
            ],
            [
                "/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/2025-12-15/electron.json.gz",
                "Electron-ID-SF",
                {
                    "year": "2024Prompt",
                    "WorkingPoint": {
                        "RecoBelow20": [0, 20],
                        "Reco20to75": [20, 75],
                        "RecoAbove75": [75, "inf"]
                    }
                }
            ]
        ],
        "muon_sf": [
            [
                "/cvmfs/cms-griddata.cern.ch/cat/metadata/MUO/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/2025-11-27/muon_Z.json.gz",
                "NUM_MediumID_DEN_TrackerMuons",
                {}
            ],
            [
                "/cvmfs/cms-griddata.cern.ch/cat/metadata/MUO/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/2025-11-27/muon_Z.json.gz",
                "NUM_LooseMiniIso_DEN_MediumID",
                {}
            ]
        ]
   }

A few practical points:

- The ID and isolation entries are independent. If you change your muon ID, you almost certainly need to change the isolation entry too,
  since the isolation SF is measured *with respect to* a specific ID (the ``DEN_TightID`` part of the correction name).
- The electron reconstruction SF uses different working-point keys for different :math:`p_T` ranges. 
  The dictionary form of ``WorkingPoint`` above selects ``RecoBelow20`` for :math:`p_T<20`\ GeV,
  ``Reco20to75`` for :math:`20<p_T<75`\ GeV, and ``RecoAbove75`` above that.
  The :doc:`/user_guide/scale_factors` page documents this range-dependent syntax in full.
- ROOT-format SF files are still common for Run 2 UL analyses. The same config keys accept either format; 
  only the inner array shape differs.


Step 4: Jet energy corrections and JME-provided selections
----------------------------------------------------------

Whether you need to reapply JEC depends on your analysis. If you trust
the JEC already applied in NanoAOD, leave ``reapply_jec`` unset and
Pepper will use the values as they are. To reapply -- which is needed
whenever the latest JEC version is more recent than the one baked into
your NanoAOD, or when you want JEC uncertainties -- set
``reapply_jec`` and configure the corrections.

For Run 3 the correctionlib-based form bundles all jet corrections
into one file:

.. code-block:: json

   {
        "reapply_jec": true,
        "smear_met": false,
        "jme_correctionlib_corrections": {
            "path": "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/2026-06-05/jet_jerc.json.gz",
            "jet_correction_data": [
                "Summer24Prompt24_V3_DATA_L1FastJet_AK4PFPuppi",
                "Summer24Prompt24_V3_DATA_L2Relative_AK4PFPuppi",
                "Summer24Prompt24_V3_DATA_L3Absolute_AK4PFPuppi",
                "Summer24Prompt24_V3_DATA_L2L3Residual_AK4PFPuppi",

            ],
            "jet_correction_mc": [
                "Summer24Prompt24_V3_MC_L1FastJet_AK4PFPuppi",
                "Summer24Prompt24_V3_MC_L2Relative_AK4PFPuppi",
                "Summer24Prompt24_V3_MC_L3Absolute_AK4PFPuppi",
            ],
            "jet_uncertainty_template": "Summer24Prompt24_V3_MC_[UNC]_AK4PFPuppi"
            "jet_uncertainty": ["Total"],
            "jet_resolution": "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi",
            "jet_ressf": "Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi",
            "jet_ressf_uncertainty": "Summer24Prompt24_JRV1_MC_SFUncertainty_AK4PFPuppi"
        }

   }

The ``[UNC]`` placeholder in ``jet_uncertainty_template`` is filled in once per entry in ``jet_uncertainty``, 
so adding more individual sources is just a matter of extending that list. 
Use ``"Total"`` for analyses that do not need the full source decomposition.

.. note::

   The JER tags above were updated in June 2026 and they contain the new ``SFUncertainty`` entry. If you are using an older JER version
   without that entry, then you should omit the ``jet_ressf_uncertainty`` key.

.. caution::

   Pepper rejects configurations that mix the legacy TXT-file keys
   (``jet_correction_mc``, ``jet_uncertainty`` at the top level, etc.)
   with ``jme_correctionlib_corrections``. Pick one approach per configuration.

Two more JME-provided files are usually configured alongside the JEC, because you will be looking them up on the same POG page anyway. 
They are not weights but event/object selections derived from POG inputs:

.. code-block:: json

   {
        "jet_ids": [
            "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/2025-12-02/jetid.json.gz",
            "AK4PUPPI"
        ],
        "jet_veto_map": [
            "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/2025-12-02/jetvetomaps.json.gz",
            "Summer24Prompt24_RunBCDEFGHI_V1",
            {"type": "jetvetomap"}
        ],
   }

``jet_ids`` set the jet ID flags that were removed in NANOAODv15; ``jet_veto_map`` removes jets falling in hot/cold detector 
regions that the POG has flagged as problematic. The veto map apply to both data and MC, and is strongly recommended for Run 3 analyses.


Step 5: B-tagging scale factors
-------------------------------

B-tagging SFs depend on per-analysis efficiencies, so they need a one-time preparation step before you can use them. 
See :ref:`bundled-scripts-example` for the full procedure -- 
the short version is that you run the main processor once with the ``btageff`` histogram enabled 
(and ``btag_sf`` *not* yet in the config), pass the resulting histogram file to ``generate_btag_efficiencies.py`` to
produce ``btag_eff.root``, and only then add the ``btag_sf`` key.

The configuration looks like:

.. code-block:: json

   {
        "btag": "upart:medium",
        "btag_wps": "$CONFDIR/btag_wps.json",
        "btag_sf": [
            [
                "/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/2026-03-10/btagging.json.gz",
                "$DATADIR/btag/btag_sf.root"
            ],
        ],
   }

The ``btag`` key selects the tagger and working point used in the selection itself; ``btag_wps`` provides the numeric thresholds; 
and ``btag_sf`` combines the central correctionlib SF file with your analysis-specific efficiency file. 
All three need to refer to the same tagger and working point -- inconsistency here produces silently wrong weights.

The ``--cut`` argument to ``generate_btag_efficiencies.py`` (which defaults to ``HasJets``) should match the cut in your selection
immediately *before* the b-tagging requirement. Mismatched cuts produce subtly wrong efficiencies that the script will not warn you about 
-- if your b-tagging weights look off, this is the first thing to check.

.. note::

   ``btag_sf`` is optional. A configuration without it runs without b-tagging weights applied (data/MC b-tag yields will not agree).

.. note::

    The ``btag_wps.json`` file is a simple dictionary mapping tagger and working point names to numeric thresholds.
    An example can be found in :repo:`example/btag_wps.json`. Pepper does not yet support reading the thresholds from correctionlib files, 
    but this is on the roadmap for a future release.

Step 6: Verifying the result
----------------------------

After each step, two quick checks catch the majority of misconfigurations:

#. **Cutflow yields.** Re-run with ``--debug`` and inspect the cutflow. 
   Total weighted yields should change roughly as expected when each correction is added 
   (a few-percent shift per lepton SF, larger for pileup, potentially larger for b-tagging in tagged regions). 
   A jump of orders of magnitude usually means a normalisation issue.
#. **Up/down systematics.** With ``"compute_systematics": true``, the output histograms should contain
   ``_up`` and ``_down`` variants for each correction you have added. 
   If a variant is missing, the corresponding correction is not actually being applied as a weight.

For a more systematic comparison against a reference set of plots, use ``make_plots.py`` (see :doc:`/user_guide/scripts`).


See also
--------

- :doc:`/user_guide/scale_factors` -- exhaustive reference for every
  SF key.
- :ref:`extending-config-example` -- if you need to add an
  analysis-specific scale factor that does not match any of the
  built-in keys.
- :ref:`bundled-scripts-example` -- the full analysis workflow,
  including the b-tagging efficiency loop.
- :repo:`example/config_ttbarll.json` -- a complete (Run 2 UL)
  configuration with all of the above corrections enabled
  simultaneously, useful as a side-by-side reference for a different
  era.