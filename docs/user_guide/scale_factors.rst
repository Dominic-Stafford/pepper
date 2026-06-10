Scale Factors
=============

On this page, you can find information about the various scale factors that Pepper supports out-of-the-box, including their definitions and usage, and any relevant corrections. 

For more detailed information, please refer to the `CMS Analysis Corrections <https://cms-analysis-corrections.docs.cern.ch/>`__ documentation
or the individual POG documentation pages.

Jet Energy Corrections and Resolution
--------------------------------------

Pepper provides two complementary pathways for applying jet energy corrections
(JEC) and jet energy resolution (JER) smearing, corresponding to the legacy
Run 2 approach using flat TXT files and the newer correctionlib-based approach
introduced for Run 3 (and hopefully available for all of Run 2 in the future).
To understand what these corrections do, please refer to the
`Jet Energy Resolution and Corrections <https://cms-jerc.web.cern.ch/>`__ documentation.

Reapplying Jet Energy Corrections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, Pepper assumes the jet corrections already applied during NanoAOD
production are sufficient. To override this and reapply corrections from
scratch, set::

   "reapply_jec": true

When this flag is set, Pepper first strips the existing NanoAOD corrections
from jets before applying the corrections you specify. This is necessary
whenever you want to use a different JEC version than the one embedded in the
NanoAOD. Setting this flag requires you to specify the corrections you want
to apply instead. See the following sections for more.

Run 2: TXT-File-Based Corrections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For Run 2 analyses, corrections are specified as arrays of paths to AK4PFchs
TXT files in the standard CMS JEC format (L1, L2, L3 levels). Commonly, 
the compounded correcions are also provided in a single txt file called 
``L1L2L3Res`` which is equivalent to the list of the individual corrections.

**Monte Carlo corrections** are given as a flat array::

   "jet_correction_mc": [
     "path/to/L1FastJet.txt",
     "path/to/L2Relative.txt",
     "path/to/L3Absolute.txt"
   ]

**Data corrections** can either be given as a flat array (if one set of
corrections covers the full dataset) or as a dictionary keyed by era name,
which allows different correction files to be applied to different run
eras::

   "jet_correction_data": {
     "RunB": ["path/to/RunB_L1.txt", "..."],
     "RunC": ["path/to/RunC_L1.txt", "..."]
   }

Jet Energy Uncertainties
"""""""""""""""""""""""""""""

Jet energy scale uncertainties are provided via a single TXT file, which can
either contain the total uncertainty or be split into individual sources::

   "jet_uncertainty": "path/to/UncertaintySources_AK4PFchs.txt"

If the file contains multiple sources, you can restrict which sources are
propagated using ``junc_sources_to_use``, supplying an explicit list of source
names. Sources not listed are ignored and all sources are included when this is unspecified.
In case you include it, your config would look like::

    "jet_uncertainty": "path/to/UncertaintySources_AK4PFchs.txt"
    "junc_sources_to_use": ["source1", "source2"]


Jet Energy Resolution
"""""""""""""""""""""""""

Legacy format (before the 05/06.2026 update)

JER smearing requires two inputs: the :math:`p_T` resolution itself and a data/MC scale
factor used to adjust the smearing applied in simulation::

   "jet_resolution": "path/to/Resolution_AK4PFchs.txt",
   "jet_ressf":      "path/to/ResolutionSF_AK4PFchs.txt"

Both must be provided together. The smearing is applied to MC jets to bring
the simulated resolution in line with what is observed in data.

New format (starting from the 05/06.2026 update)

SF and SF uncertainty are provided as distinct inputs. So the new JER smearing configuration looks like::

    "jet_resolution":"Summer24Prompt24_JRV1_MC_PtResolution_AK4PFPuppi",
    "jet_ressf":"Summer24Prompt24_JRV1_MC_ScaleFactor_AK4PFPuppi",
    "jet_ressf_uncertainty": "Summer24Prompt24_JRV1_MC_SFUncertainty_AK4PFPuppi"

MET Smearing
"""""""""""""""""""""""""""

When jet momenta are modified by smearing, the missing transverse energy (MET) can
be updated accordingly. This propagation is controlled by::

   "smear_met": true

Run 3: Correctionlib-Based Corrections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Starting from Run 3, the preferred approach consolidates all JME corrections
into a single correctionlib JSON file. This is configured via
``jme_correctionlib_corrections``, which takes a dictionary with the following
structure::

   "jme_correctionlib_corrections": {
     "path": "path/to/corrections.json.gz",
     "jet_correction_data": [
                                "L1FastJet_DATA", 
                                "L2Relative_DATA", 
                                "L3Absolute_DATA"
                            ],
     "jet_correction_mc": [
                            "L1FastJet_MC",
                            "L2Relative_MC",
                            "L3Absolute_MC"
                          ],
     "jet_uncertainty_template": "Campaign_V3_MC_[UNC]_AK4PFPuppi",
     "jet_uncertainty": [
                            "AbsoluteStat",
                            "AbsoluteScale",
                            "FlavorQCD"
                        ],
     "jet_resolution": "JetResolution_AK4PFPuppi",
     "jet_ressf": "JetResolutionSF_AK4PFPuppi"
   }

The ``jet_correction_data`` and ``jet_correction_mc`` keys list the named
correction levels within the correctionlib file that should be compounded and
applied in order. 
The uncertainty template uses the ``[UNC]`` wildcard, which
Pepper substitutes with each entry in ``jet_uncertainty`` to locate the
corresponding uncertainty object inside the file. This makes it straightforward
to add or remove uncertainty sources without changing the template string.
Commonly, an uncertainty named ``Total`` is included in the correctionlib file
that encompasses all sources.  Commonly, the compounded correcions are also 
provided in a single txt file called ``L1L2L3Res`` which is equivalent to the 
list of the individual corrections.

Summary of Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The table below summarises which inputs are required or implied by each
configuration choice.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Goal
     - Required keys
   * - Use NanoAOD JEC as-is
     - *(nothing)*
   * - Reapply JEC (Run 2)
     - ``reapply_jec``, ``jet_correction_mc``, ``jet_correction_data``
   * - JEC uncertainties (Run 2)
     - ``jet_uncertainty`` (requires ``jet_correction_mc``)
   * - JER smearing (Run 2)
     - ``jet_resolution``, ``jet_ressf``
   * - All of the above (Run 3)
     - ``reapply_jec``, ``jme_correctionlib_corrections``
   * - Propagate smearing to MET
     - ``smear_met: true``

.. note::

   ``jet_correction_mc`` is also required whenever ``jet_uncertainty`` or
   ``jet_resolution`` is specified, even if ``reapply_jec`` is false, since
   the correction levels are needed to define the reference point for
   uncertainty and smearing calculations.


Jet Veto Maps
--------------

Jet veto maps allow specific regions of the detector to be masked out by
vetoing jets that fall within known problematic areas
(see more `here <https://cms-jerc.web.cern.ch/Recommendations/#jet-veto-maps>`__).
They are configured via ``jet_veto_map``, an optional array of arrays.
Each inner array contains two elements:

1. A path to a correctionlib JSON file containing the veto map.
2. A two-element array specifying the campaign/era name and the map name within
   the file.

For example::

   "jet_veto_map": [
     [
       "path/to/jetvetomap.json.gz",
       ["Campaign_Name_RunEFG_V1", "jetvetomap"]
     ]
   ]

Currently, Pepper only supports a single veto map in the list.

Pile-Up
----------
Pileup reweighting corrects for differences in the pileup distribution between
data and simulation. It is configured via ``pileup_reweighting``, which accepts
either a ROOT file or a correctionlib-based list.

**Correctionlib list** (preferred)

Pass a list defining a correctionlib correction. The elements to be specified are:

1. Path to the correctionlib JSON file.
2. The correction key identifying the specific pileup profile to use (e.g. ``"Collisions17_UltraLegacy_goldenJSON"``).
3. *(Optional)* A dictionary of extra parameters required by the correctionlib recipe (typically omitted for pile-up corrections).

For example::

   "pileup_reweighting": [
       "path/to/puWeights.json.gz",
       "Collisions17_UltraLegacy_goldenJSON"
     ]

**ROOT file** (Legacy approach)

Pass the path to a ROOT file produced by ``compute_pileup_weights.py``::

   "pileup_reweighting": "path/to/pileup_weights.root"


Photon Scale Factors
---------------------

Photon Scale Smearing
---------------------

Electron ID and Reco Scale Factors
------------------------------------

The electron scale factors - reconstruction (Reco), identification (ID), and
isolation (ISO) - are configured through the single ``electron_sf`` key. Its
value is an optional array of arrays, where each inner array represents one
correction. The total electron weight applied to an event is the product of all
entries, so stacking Reco, ID, and ISO corrections is simply a matter of
listing them in order.

Each inner array contains exactly three elements and can refer to either a
ROOT-based or correctionlib-based source.

**ROOT format**

Provide the path to a ROOT file, the name of the 2D histogram containing the
scale factors, and the axis order as a two-element array::

   "electron_sf": [
     ["path/to/egamma_reco_sf.root", "EGamma_SF2D", ["eta", "pt"]],
     ["path/to/egamma_id_sf.root",   "EGamma_SF2D", ["eta", "pt"]]
   ]

The axis order must match the binning of the histogram. The typical histogram
name is ``EGamma_SF2D``.

**Correctionlib format**

Provide the path to a correctionlib JSON file, the correction name, and a
dictionary of any additional input variables required by the correction (such
as the year or working point)::

   "electron_sf": [
     [
       "path/to/electron_sf.json.gz",
       "UL-Electron-ID-SF",
       {"year": "2017", "WorkingPoint": "RecoAbove20"}
     ],
     [
       "path/to/electron_sf.json.gz",
       "UL-Electron-ID-SF",
       {"year": "2017", "WorkingPoint": "wp80iso"}
     ]
   ]

**pT-range-dependent correction keys**

Some corrections, such as the reconstruction SF, use different working point
keys for different :math:`p_T` ranges as documented in the 
`EGamma POG <https://twiki.cern.ch/twiki/bin/view/CMS/EgammSFandSSRun3>`__.
In this case, the working point value in the extra-info dictionary can itself 
be a dictionary mapping key names to ``[lo, hi]`` :math:`p_T` intervals.
The upper bound of the last interval should be ``"inf"``::

   "electron_sf": [
     [
       "path/to/electron_sf.json.gz",
       "UL-Electron-ID-SF",
       {
         "year": "2017",
         "WorkingPoint": {
           "Reco20to75":  [20, 75],
           "RecoAbove75": [75, "inf"]
         }
       }
     ],
     [
       "path/to/electron_sf.json.gz",
       "UL-Electron-ID-SF",
       {"year": "2017", "WorkingPoint": "wp80iso"}
     ]
   ]

This :math:`p_T`-splitting mechanism is particularly common for the Reco SF, which the
EGamma POG provides in separate low-:math:`p_T` and high-:math:`p_T` variants.


Electron Scale Smearing
-------------------------

Electron Trigger Scale Factors
----------------------------------

Muon ID and ISO Scale Factors
-------------------------------

Muon identification (ID) and isolation (ISO) scale factors are both configured
through the single ``muon_sf`` key, following the same structure as
``electron_sf``. Each entry in the outer array represents one correction, and
the total muon weight is the product of all entries. ID and ISO scale factors
are therefore applied by listing them as separate entries.

Each inner array contains exactly three elements and supports both ROOT-based
and correctionlib-based sources.

To find more information about the available muon scale factors, consult the
`Muon POG <https://muon-wiki.docs.cern.ch/>`__.

**ROOT format**

Provide the path to a ROOT file, the histogram name, and the axis order. For
muons the typical axis order is ``["abseta", "pt"]``, reflecting that muon SFs
are usually binned in absolute pseudorapidity::

   "muon_sf": [
     ["path/to/muon_id_sf.root",  "NUM_MediumID_DEN_genTracks_abseta_pt",  ["abseta", "pt"]],
     ["path/to/muon_iso_sf.root", "NUM_TightRelIso_DEN_MediumID_abseta_pt", ["abseta", "pt"]]
   ]

**Correctionlib format**

Provide the path to a correctionlib JSON file, the correction name, and a
dictionary of any additional input variables required by the correction::

   "muon_sf": [
     [
       "path/to/muon_sf.json.gz",
       "NUM_MediumPromptID_DEN_TrackerMuons",
       {}
     ],
     [
       "path/to/muon_sf.json.gz",
       "NUM_TightRelIso_DEN_MediumPromptID",
       {}
     ]
   ]

**Splitting systematic and statistical uncertainties**

Muon SF uncertainties can optionally be decomposed into a systematic component
(correlated across years) and a statistical component (uncorrelated across
years) by setting::

   "split_muon_uncertainty": true

This splits the muon uncertainty into statistical and systematic components, as provided by the Muon POG.


Muon Scale Smearing
--------------------

Muon Rochester
--------------------

Tau Scale Factors
------------------

Flavor Tagging Scale Factors
-----------------------------

B-tagging scale factors correct for differences in tagging efficiency between data and simulation.

The b-tagging corrections are specified under the key ``btag_sf`` and require a ROOT file calculated on a
per-analysis level using the script ``generate_btag_efficiencies.py``.
Refer to `BTV POG documentation <https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods>`__ for more
details on how b-tagging efficiency corrections work.

The resulting configuration entry looks like::

  "btag_sf": [
      [
          "path/to/btagging.json.gz",
          "path/to/btag_sf.root"
      ],
  ],

To calculate ``btag_sf.root``, use the script ``generate_btag_efficiencies.py``.
The script expects a ``hists.json`` file produced by your analysis processor and a ``--cut`` argument
specifying the name of the cut immediately before the b-tagging requirements (default ``'HasJets'``).
