.. _configuration-reference:

Pepper Configuration
====================

Pepper uses a JSON-based configuration format. JSON arrays map to Python
lists, and JSON objects map to Python dictionaries with string keys.
A variable set to ``null`` in the JSON file is treated as if it were not
present.

In addition to strict JSON, Pepper also accepts **HJSON** as an input format,
which allows comments, multi-line strings, unquoted keys, and other
user-friendly features. See the HJSON documentation for details:

   https://hjson.github.io/

If a configuration variable expects a *path to another JSON/HJSON file*, its
value may also be **the file’s JSON/HJSON content itself**, inlined directly.

Masses and transverse momenta (``pT``) are always interpreted in **GeV**.

The following sections describe the available configuration options. Each of the 
options corresponds to a key in the JSON configuration file. Options marked as
*optional* may be omitted.

--------------------
Special Variables
--------------------

``import``
   Optional. Path to another configuration JSON/HJSON file for inheritance. That file is
   loaded *before* the current one, and values in the current file override
   values from the imported file.

The following string placeholders may appear inside configuration values and
are substituted at runtime:

``$CONFDIR``
   Path to the directory containing the configuration file supplied on the
   command line.

.. _storedir:

``$STOREDIR``
   The value of the :ref:`store <config-store>` configuration key.

``$DATADIR``
   The value of the :ref:`datadir <config-datadir>` configuration key.

--------------------
General Options
--------------------

``year`` (string)
   The data-taking year. Some settings (e.g., b-tagging working points) depend
   on this value.

``sqrt_s`` (float, optional)
   Center-of-mass energy in GeV.

``rng_seed_file`` (path, optional)
   Path to a text file used to store and load an integer seed for the random
   number generator. If the file exists, the stored seed is used; otherwise a
   seed is created and written. The final seed is combined with a
   chunk-specific offset so that each event chunk receives unique randomness.

``blinding_denom`` (float, optional)
   Blinding factor for data. Only ``1 / blinding_denom`` of data events is used.
   Monte Carlo event weights are scaled accordingly.

``compute_systematics`` (boolean)
   If ``true``, compute all systematic variations.

``skip_nonshape_systematics`` (boolean, optional)
   If ``false`` and ``compute_systematics`` is ``true``, histograms are also
   produced for non-shape uncertainties such as luminosity. If ``true``, only
   shape uncertainties produce histograms.

``exit_on_failed_jobs`` (string, optional)
   Policy for handling jobs that fail all retry attempts. Allowed values:

   - ``"all"`` – Exit on any job failure (recommended when using
     pre-computed lumifactors).
   - ``"data"`` – Exit only on data job failures (useful with
     post-computed lumifactors).
   - ``"none"`` – Never exit automatically; the user is responsible for
     handling failures.


-------------
Data Sets
-------------

``file_mode`` (string, optional)
   Determines how dataset names are resolved. Allowed values:

   - ``"local"`` (default) – Files from datasets are searched for and opened
     only from subdirectories of ``store``.
   - ``"xrootd"`` – Datasets are resolved online and files are accessed
     via XRootD.
   - ``"local+xrootd"`` – Files are first searched for in ``store``; if not
     found, they are accessed via XRootD.

   For direct paths or glob patterns, this option has no effect.

``xrootddomain`` (string)
   Domain (redirector) used for XRootD requests.

``file_blacklist`` (array or path, optional)
   A list, or a path to a JSON/HJSON file containing a list, of files or
   logical file names to skip entirely.

``local_file_blacklist`` (array or path, optional)
   A list, or a path to a JSON/HJSON file, containing paths or logical file
   names that should always be resolved remotely via XRootD, even if they
   exist locally. Only applicable when ``file_mode = "local+xrootd"``.

``xrootd_url_blacklist`` (array, optional)
   List of XRootD URLs or URL fragments (e.g., ``"rl.ac.uk"`` or ``".fr"``)
   that should not be used when resolving files via XRootD.

.. _config-store:

``store`` (string)
   Base directory used to resolve dataset names and substituted into
   ``$STOREDIR``.

``exp_datasets`` (object)
   Maps dataset labels to arrays of dataset sources. Each array contains:

   - Paths or XRootD URLs to NanoAOD ROOT files (wildcards allowed), or  
   - Dataset names in the CMS three-slash format (``/A/B/C``).

   These datasets are used for experimental data input.

``MET_trigger_datasets`` (object)
   Same structure as ``exp_datasets``. Used as cross-trigger inputs for
   trigger scale factor computation.

``data_eras`` (object)
   Maps era names to arrays containing two-element lists specifying the run
   number ranges corresponding to each era.

``mc_datasets`` (object)
   Same structure as ``exp_datasets``. Provides Monte Carlo dataset inputs.

``dataset_for_systematics`` (object)
   Maps MC dataset keys to arrays of two strings:

   1. The name of the nominal MC dataset (a key in ``mc_datasets``).  
   2. The name of the systematic uncertainty and its direction.

   Indicates that the dataset represents a systematic variation of another.

``dataset_trigger_map`` (object)
   Maps dataset names (keys from ``exp_datasets``) to arrays of trigger
   paths. If a trigger is absent from the NanoAOD file, it is ignored.  
   The user must ensure there are no typos.

``MET_triggers`` (array of strings)
   Trigger paths used as cross-triggers for MET trigger scale factors.

``dataset_trigger_order`` (array)
   Lists keys from ``exp_datasets`` in priority order. If an event passes
   multiple dataset trigger selections, it is processed only in the
   first dataset listed here.

``crosssection_uncertainty`` (path or object, optional)
   Path to a JSON/HJSON file containing an object of arrays, or the object
   itself. Keys correspond to MC dataset names. Each array contains:

   1. A string naming the cross-section uncertainty.
   2. A float giving the uncertainty factor, or ``null`` to indicate no
      cross-section uncertainty.

   Two uncertainties with the same name are treated as fully correlated.
   MC datasets not listed here are assigned zero cross-section uncertainty.
   This does *not* affect datasets listed in ``dataset_for_systematics``.

``pdf_types`` (object, optional)
   Maps PDF set LHA IDs to their type. Allowed values:

   - ``"Hessian"``
   - ``"MC_Gaussian"`` – Use standard deviation of replicas.
   - ``"MC"`` – Use ordered distribution and extract the 16%/84% envelope.

   Pepper’s example configuration covers most Run 2 datasets; more can be
   added as needed.

``split_pdf_uncs`` (bool, optional)
   Output the full set of PDF variations instead of computing the combined
   uncertainty. For NLO samples, this is usually required for accurate
   histogram-level uncertainty evaluation.

``normalize_pdf_uncs`` (bool, optional)
   Normalize PDF variations to the nominal cross section before cuts. This
   removes cross-section-only effects while preserving shape changes and
   acceptance-induced rate effects. Requires running
   ``compute_mc_lumifactors.py`` with ``-p``.

``DY_datasets`` (array, optional)
   List of MC datasets used for DY scale factor computation and application.
   If omitted, all MC datasets starting with ``"DY"`` are used.

------------------------------------------------
Scale Factors, Weights, and Other Inputs
------------------------------------------------

These options determine various calibrations and weightings. If an
optional parameter is missing, the corresponding procedure will be skipped.

.. _config-datadir:

``datadir`` (string, optional)
   Base directory used for data inputs. Substituted into ``$DATADIR``.

``lumimask`` (path, optional)
   Path to the Golden JSON file specifying certified good runs.

``luminosity`` (float)
   Integrated luminosity in 1/fb. Required when computing
   ``mc_lumifactors``. Represents the luminosity after applying the
   ``lumimask``.

``crosssections`` (object or path)
   Object containing floats or a path to a JSON/HJSON file with floats.
   Required for ``mc_lumifactors`` computation. Provides cross sections
   for all datasets in ``mc_datasets`` (units: fb).

``norm_genweights`` (bool, optional)
   Whether to normalize all generator weights, considering only the sign.
   Recommended if a sample has a few very large event weights.

``genweights_to_norm`` (array, optional)
   If ``norm_genweights`` is false, normalize only this list of datasets.

``mc_lumifactors`` (path or false)
   Path to a JSON file produced by ``compute_mc_lumifactors.py`` or
   ``false``. If ``false``, Pepper will compute lumifactors at runtime.

``pileup_reweighting`` (path or array (of length 2 or 3), optional)
   ROOT file containing pileup weights produced by
   ``compute_pileup_weights.py`` or a list of correctionlib corrections. For each item in a correctionlib list,
   the first element is the path to the json, the second is the key for the corrections (e.g. ``"Collisions17_UltraLegacy_goldenJSON"``),
   and the third, which can be ommited, is a dictionary for potential extra parameters required by the corrections.

``stitching_factors`` (object, optional)
   Dictionary mapping dataset roots to stitching factors. Each factor entry
   should include:

   - ``axis`` – Datapicker object specifying the axis for stitching.  
   - ``edges`` – Array defining bin edges (last edge = infinity).  
   - ``factors`` – Output from ``calculate_stitching_factors.py``.

``electron_sf`` (array of arrays, optional)
   Each inner array must contain exactly three elements. Can refer to ROOT
   SFs or correctionlib SFs:

   - ROOT: [path_to_ROOT, histogram_name, axis_order].  
     Typical histogram_name: ``EGamma_SF2D``.  
     Axis order: e.g., ``["eta", "pt"]``.
   - correctionlib: [path_to_json, correction_name, extra_info].  
     Example: ``"year"`` and ``WorkingPoint``.
   - Can specify different correction keys for :math:`p_T` ranges using a dictionary:
     ``"WorkingPoint": {"Reco20to75": [20, 75], "RecoAbove75": [75, "inf"]}``.

``muon_sf`` (array of arrays, optional)
   Same structure as ``electron_sf``. ROOT axes typically ``["abseta", "pt"]``.

``split_muon_uncertainty`` (bool, optional)
   If true, split muon SF uncertainties into systematic and statistical parts.
   Systematic uncertainties can be treated as correlated across years.

``muon_rochester`` (path, optional)
   Path to Muon POG Rochester correction TXT file.

``btag_sf`` (array of arrays, optional)
   Each inner array contains exactly two paths:

   1. CSV file with BTV POG scale factors.  
   2. ROOT file with efficiencies from ``generate_btag_efficiencies.py``.

``btag_method`` (string, optional)
   Method for b-tagging SF computation: ``fixedwp`` or ``iterativefit``.
   Default: ``fixedwp``.

``btag_measure_type`` (string, optional)
   SF measurement type: ``mujets`` (muon-enriched QCD) or ``comb`` (combined
   samples, e.g., ttbar). Default: ``mujets``.

``btag_ignoremissing`` (bool, optional)
   If true, missing SFs return 1 instead of throwing an exception. Default: false.

``btag_splitting_scheme`` (string, optional)
   Splitting scheme for b-tag uncertainties:

   - ``Null`` – no splitting  
   - ``years`` – correlated and uncorrelated across years  
   - ``sources`` – split into all available subsources

``split_btag_year_corr`` (bool, optional, deprecated)
   Replaced by ``btag_splitting_scheme``. Splits btag uncertainties across years.

``jet_veto_maps`` (array of arrays, optional)
   Each inner array must contain a string and a two-tuple. The string is the path to the JSON file with jet veto maps, and the two-tuple contains the era and the map name::

      "jet_veto_maps": [
         [
               "path/to/jetvetomaps.json.gz",
               ["Campaign_Name_RunEFG_V1", "jetvetomap"],
         ]
      ],

``jet_puid_sf`` (array, optional)
   First element: path to Jet PU ID SF JSON (correctionlib).  
   Second element (optional): path to ROOT MC efficiencies.

``reapply_jec`` (bool, optional)
   If true, undo existing NanoAOD jet corrections and reapply
   ``jet_correction`` and ``jet_correction_data``.

``jme_correctionlib_corrections`` (dict, optional)
   Starting in Run 3, this is the prefered way to define JME corrections using a single correctionlib file.
   Needed if ``reapply_jec`` is true. Must point to a dictionary with the following keys:

   - ``path``: Path to the correctionlib json file
   - ``jet_correction_data``: Array of corrections within file to compound and apply to data
   - ``jet_correction_mc``: Array of corrections within file to compound and apply to MC
   - ``jet_uncertainty_template``: A template string for how the uncertainties are called. In Run 3 2023Bpix this could be ``Summer23BPixPrompt23_V3_MC_[UNC]_AK4PFPuppi``. Notice the wildcard ``[UNC]`` which will be substituted with the items in ``jet_uncertainty`` to find the specific uncertainty.
   - ``jet_uncertainty``: Array of strings. Each string is the name of a jet uncertainty to be applied. The name is substituted into ``jet_uncertainty_template``.
   - ``jet_resolution``: A string with the jet resolution correcion name.
   - ``jet_ressf``: A string with the jet resolution scale factor name.

``jet_correction_mc`` (array)
   Array of paths to AK4PFchs MC JEC TXT files (L1, L2, L3). Needed if
   ``jet_uncertainty`` or ``jet_resolution`` is given, or if ``reapply_jec`` is true.

``jet_correction_data`` (array or dict)
   Paths to AK4PFchs data JEC TXT files. Can be split by era with a dict.

``jet_uncertainty`` (path, optional)
   Path to AK4PFchs jet energy uncertainty TXT file. Can be total or per-source.

``junc_sources_to_use`` (array of strings, optional)
   List of uncertainty sources to use instead of all sources in
   ``jet_uncertainty``.

``jet_resolution`` (path, optional)
   AK4PFchs jet pt resolution TXT file (used with ``jet_ressf``).

``jet_ressf`` (path, optional)
   Jet resolution SF TXT file, used with ``jet_resolution``.

``smear_met`` (bool, optional)
   Whether to propagate jet smearing to MET.

``MET_xy_shifts`` (path, optional)
   Path to JSON file from ``produce_met_xy_nums.py``.

``top_pt_reweighting`` (object, optional)
   Defines top pT reweighting. Recognized keys:

   - ``method`` – ``"theory"`` or ``"datanlo"``  
   - For ``"theory"``: ``a, b, c, d`` (float)  
     Formula: ``sf = a * exp(b*pt) - c*pt + d``  
   - For ``"datanlo"``: ``a, b`` (float)  
     Formula: ``sf = exp(a + b*pt)``  
   - ``scale`` (float, optional) – global factor  
   - ``sys_only`` (bool, optional) – apply as systematic only, without
     multiplying nominal weight

--------------------
Output
--------------------

``columns_to_save`` (array or object, optional)
   Data pickers specifying which observables to save per event.  
   If an object, the keys are used inside the output file for the corresponding observable.  
   Note: Per-event output must be explicitly enabled even if this variable is present.

``column_output_format`` (string, optional)
   File format used to save the columns in ``columns_to_save``.  
   Options: ``"root"``, ``"hdf5"``. Default: ``"root"``.

``save_categories_per_event`` (bool, optional)
   Whether to save all categories in per-event output. Default: ``True``.  
   If ``False``, individual columns must be specified in ``columns_to_save``.

``use_temp_eventdir`` (bool, optional)
   Write out the per-event output in a temporary directory on the worker, and then copy it to the final destination
   afterwards. On CERN EOS, ``eos cp`` is used for the copying instead of the file system mount, 
   which should be more stable. Default is ``true`` on EOS and ``false`` otherwise.

``hists`` (path or object)
   Path to a JSON/HJSON file containing histogram definitions as an object.  
   Keys determine histogram names. See *Histogram Definition* below.

``cuts_to_histogram`` (array of strings, optional)
   Only create histograms for cuts whose names are in this array.

``systs_to_histogram`` (array of strings, optional)
   Only include systematics listed here. For two-sided or numeric systematics,
   do **not** include direction/indices.  
   Example: use ``"muonsf"``, not ``"muonsf_down"``.

``datasets_to_group`` (object, optional)
   Maps dataset names to group names for histograms.

``histogram_format`` (string, optional)
   Output histogram format: ``"hist"``, ``"coffea"``, or ``"root"``. Default: ``"hist"``.

^^^^^^^^^^^^^^^^^^^^
Histogram Definition
^^^^^^^^^^^^^^^^^^^^

Elements of the ``hists`` object are objects themselves, with the following keys:

``bins`` (array of objects, optional)
   Defines binned axes. Each axis object contains:

   - ``name`` – Name of the axis  
   - ``label`` – Axis label (supports LaTeX with dollar signs)  
   - ``n_or_arr`` – Number of bins (regular) or array of bin edges  
   - ``lo`` – Minimum edge (if regular)  
   - ``hi`` – Maximum edge (if regular)  
   - ``unit`` – Optional string specifying unit (e.g., ``"GeV"``)  
   - ``type`` – Optional, ``float`` or ``int`` (default: ``float``)

``cats`` (array of objects, optional)
   Category axes to group events. Each object contains:

   - ``name`` – Name of the category axis  
   - ``label`` – Label for plotting  
   - ``default_key`` (optional) – Default value for events missing a category

Note: Category axes for dataset and channel are automatically created.

``fill`` (object)
   Defines what to fill in the histogram. Keys must be either:

   - A key in ``bins`` → value is a data picker  
   - A key in ``cats`` → value is an object mapping category names to data pickers

If ``fill`` is missing, the histogram will not be filled.

``step_requirement`` (string, optional)
   Histogram is only created after a specific step. Format:

   - ``"cut:<cut_name>"``  
   - ``"column:<column_name>"``

``do_systs`` (bool, optional)
   Include systematics for this histogram. Default: ``true``.

``weight`` (data picker, optional)
   Custom event weight. If absent, uses default event weight.

``weight_factor`` (data picker, optional)
   Cannot be used in conjuction with ``weight``. Specifies a custom event weight 
   factor using a data picker, which will multiply the event weight (including systematics).
   Can either be one value per event, or have the same shape as the fill.

``label`` (string, optional)
   Label for the bin-height axis. Default chosen to match CMS guidelines
   (e.g., ``"Events / bin"``).

``selector_cats`` (bool or array, optional)
   Include categories from the selector axis:

   - ``True`` → include all selector categories
   - List → include only specified categories. 
   
   Default: ``true``.

``only_for_cats`` (dict, optional)
   A dictionary of the form ``{"category_name": ["category_value", ...]}``. If given,
   fill the histogram only for the given categories (e.g. to save memory). For example,
   giving ``{"channel": ["is_ee", "is_mm"]}`` would fill the histogram only for the ``ee`` and ``mumu``
   channels. If only a single value is given, encasing it in a list is optional (i.e. ``{"channel": "is_em"}`` is accepted).

^^^^^^^^^^^^^^^^^^^^
Data Pickers
^^^^^^^^^^^^^^^^^^^^

Data pickers are arrays specifying how to select data from the processor's
``data`` array. Each element refines the previous element’s output. Elements can be:

- Column name or object: ``{"key": columnname}``  
- Attribute name or object: ``{"attribute": attribute}``  
- Method name (executes and uses return value)  
- Function object: ``{"function": functionname}`` (functions defined in ``hist_defns.py``)  
- Leading selector: ``{"leading": integer}`` or ``{"leading": [integer1, integer2]}``  
  Equivalent to Python slicing: ``[:, integer - 1]`` or ``[:, integer1-1:integer2-1]``  

------------------------------
Particle Definitions and Cuts
------------------------------

Lepton and Jet selection options. These control which particles are considered
"good" for analysis, including kinematic cuts, ID requirements, isolation,
and b-tagging.

^^^^^^^^^^^^^^^^^^^^^^^^^^
MET and Electron Filters
^^^^^^^^^^^^^^^^^^^^^^^^^^

``apply_met_filters`` (bool)
   If ``true``, events are required to pass MET filters.

``ele_cut_transreg`` (bool)
   If ``true``, electrons are required to be outside the barrel-endcap
   transition region.

``ele_eta_min`` (float)
   Minimum pseudorapidity for electrons.

``ele_eta_max`` (float)
   Maximum pseudorapidity for electrons.

``good_ele_id`` (string)
   Electron ID requirement. Options:

   - ``"cut:<WP>"`` – cut-based working point, ``loose``, ``medium``, ``tight``  
   - ``"mva:<WP>"`` – MVA-based working point, ``Iso80``, ``Iso90``  
   - ``"skip"`` – no ID required

``good_ele_pt_min`` (float)
   Minimum transverse momentum for electrons.

``additional_ele_id`` / ``additional_ele_pt_min``
   Same as ``good_ele_id`` / ``good_ele_pt_min``, used for identifying
   additional leptons beyond the primary ones.

^^^^^^^^^^^^^^^^^^^^^^^^^^
Muon Selection
^^^^^^^^^^^^^^^^^^^^^^^^^^

``muon_cut_transreg`` (bool)
   If ``true``, muons are required to be outside the barrel-endcap transition
   region.

``muon_eta_min`` / ``muon_eta_max`` (float)
   Minimum and maximum pseudorapidity for muons.

``good_muon_id`` (string)
   Muon ID requirement. Can be ``cut:`` or ``mva:`` with corresponding working points.

``good_muon_iso`` (string)
   Muon isolation requirement. Options:

   - ``cut:<WP>`` – ``very_loose``, ``loose``, ``medium``, ``tight``, 
     ``very_tight``, ``very_very_tight``
   - Custom isolation: ``dR<0.3_chg:<val>``, ``dR<0.3_all:<val>``, 
     ``dR<0.4_all:<val>``

``good_muon_pt_min`` (float)
   Minimum transverse momentum for muons.

``additional_muon_id`` / ``additional_muon_iso`` / ``additional_muon_pt_min``
   Same as primary muon cuts, for identifying additional leptons. Use
   ``"skip"`` to disable ID/isolation.

^^^^^^^^^^^^^^^^^^^^^^^^^^
Jet Selection
^^^^^^^^^^^^^^^^^^^^^^^^^^

``good_jet_id`` (string)
   Jet ID requirement. ``cut:<WP>`` or ``"skip"``.

``good_jet_puId`` (string)
   Jet pileup ID requirement. ``cut:<WP>`` or ``"skip"``.

``good_jet_lepton_distance`` (float)
   Minimum ΔR distance between jets and any lepton.

``good_jet_eta_min`` / ``good_jet_eta_max`` (float)
   Minimum and maximum jet pseudorapidity.

``good_jet_pt_min`` (float)
   Minimum jet transverse momentum.

``hem_cut_if_ele`` / ``hem_cut_if_muon`` / ``hem_cut_if_jet`` (bool)
   Ignore objects in the HEM 2018 problematic region.

^^^^^^^^^^^^^^^^^^^^^^^^^^
Event-Level Cuts
^^^^^^^^^^^^^^^^^^^^^^^^^^

``mll_min`` (float)
   Minimum invariant mass for the dilepton system.

``lep_pt_min`` (array of floats)
   Minimum pT for the n-th leading lepton.

``lep_pt_num_satisfied`` (int)
   Minimum number of leptons that must satisfy the ``lep_pt_min`` cuts.

``num_jets_atleast`` (int)
   Minimum number of jets required.

``jet_pt_min`` (array of floats)
   Minimum pT for n-th leading jets.

``jet_pt_num_satisfied`` (int)
   Minimum number of jets satisfying the ``jet_pt_min`` requirements.

^^^^^^^^^^^^^^^^^^^^^^^^^^
b-Tagging
^^^^^^^^^^^^^^^^^^^^^^^^^^

``btag`` (string)
   Algorithm and working point. Options: ``deepcsv:<WP>``, ``deepjet:<WP>``.

``num_atleast_btagged`` (int)
   Minimum number of b-tagged jets required.

``btag_wp`` (object or path, optional)
   Defines b-tag working points. Can be an object or JSON file.
   Should include algorithm names, years, WP names, and one-letter abbreviations.  
   Example file: ``example/btag_wps.json``.

--------------------------
Kinematic Reconstruction
--------------------------

``reco_info_file`` (path)
   Required if ``reco_algorithm`` is present. Path to a ROOT file produced by
   ``compute_kinreco_hists.py``.

``reco_w_mass`` (float or string)
   Required if ``reco_algorithm`` is ``"Sonnenschein"``. W boson mass used
   in Sonnenschein reconstruction. If string (e.g., ``"mw"``), the value is
   randomly drawn from a histogram inside ``reco_info_file`` named by the string.

``reco_t_mass`` (float or string)
   Required if ``reco_algorithm`` is ``"Sonnenschein"``. Top quark mass used
   in Sonnenschein reconstruction. Works analogously to ``reco_w_mass``.

``reco_num_smear`` (int)
   Required if ``reco_algorithm`` is ``"Sonnenschein"``. Number of smearings
   performed in Sonnenschein.

^^^^^^^^^^^^^^^^^^^^^^^
ttbarll Specific
^^^^^^^^^^^^^^^^^^^^^^^

``channel_trigger_map`` (object of arrays)
   Keys are channels: ``"ee"``, ``"emu"``, ``"mumu"``, ``"e"``, ``"mu"``.  
   Arrays list trigger paths. Events not passing a listed trigger are ignored.

``drellyan_sf`` (path, optional)
   Drell-Yan scale factors. Can be:

   - ROOT file from ``calculate_DY_SFs.py``  
   - JSON file with object containing keys: ``bins``, ``factors``, ``factors_up``, ``factors_down``

``trigger_sfs`` (array, optional)
   First element: ROOT file containing trigger SF histograms.  
   Second element: array of three histogram names for ``ee``, ``eµ``, and ``µµ`` channels.  
   Histograms binned by leading and subleading lepton pT.

``reco_algorithm`` (string, optional)
   Algorithm for top quark reconstruction. Options: ``"Sonnenschein"``, ``"Betchart"``.

``z_boson_window_start`` / ``z_boson_window_end`` (float)
   Lower and upper bounds of invariant mass used for Z window cut.

``ee/mm_min_met`` (float)
   Minimum MET pT in the corresponding channel cut.

-----------------------------
Drell-Yan Scale Factors
-----------------------------

``fast_dy_sfs`` (bool, optional)
   If ``true``, DYprocessor runs over only DY and observed data. Not relevant for other processors.

``bin_dy_sfs`` (data picker, optional)
   Variable used to bin DY SFs when applying them via the Processor.

-----------------------------
Plotting
-----------------------------

General-purpose plotting script: ``scripts/make_plots.py``.  
Designed to be flexible and extendible for analysis-specific needs.

Plotting is controlled via a separate config file. Example: ``example/example_plotting.json``

``backgrounds`` (dict)
   Defines background processes to plot in stack plots. Format:

   - Key: process name  
   - Value: object with keys:

     - ``datasets`` – list of dataset names  
     - ``label`` – optional legend label (supports LaTeX syntax)  
     - ``color`` – optional color (name or hex code). Default: CMS official scheme

``data`` (list, optional)
   Observed datasets to plot as black markers. Default: only MC is plotted.

``signals`` (dict, optional)
   Signal processes to plot as lines, same format as ``backgrounds``.

``stack_signals`` (bool, optional)
   If ``true``, plot sum of signals and background; else plot signals separately. Default: false.

``sum_categories`` (bool, optional)
   If ``true``, sum over all categorical axes, creating one plot per histogram. Default: false.

``categories`` (dict, optional)
   Custom categories to plot. Format:

   ``"category_name": {"ax1": ["val1","val2"], "ax2": "val3", ...}``  
   - Axes not specified are summed over  
   - Ignored if ``sum_categories`` is true  
   - Default: plot all combinations of categorical axes

``rebin`` (dict, optional)
   Define new binning for specific histograms:

   - ``n_or_arr`` – integer number of bins or array of edges (required)  
   - ``lo`` – lower edge (if ``n_or_arr`` is int)  
   - ``hi`` – upper edge (if ``n_or_arr`` is int)  
   - ``label`` – optional x-axis label  

   See ``example/example_plotting.json`` for examples.

``year`` (int, optional)
   Year to display on the top right of the plot.

``lumi`` (float, optional)
   Luminosity to display.

``com_energy`` (float, optional)
   Center-of-mass energy in TeV. Default: 13
