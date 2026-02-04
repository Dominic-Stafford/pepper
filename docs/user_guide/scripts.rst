Main Scripts
============

The ``scripts`` directory contains helper scripts to obtain inputs and produce or visualize outputs:

- ``calulate_stitching_factors.py``  
  Calculate stitching factors from an already produced histogram  
  (currently only 1D histograms supported).

- ``compute_mc_lumifactors.py``  
  Compute :math:`\mathcal{L}\sigma / \sum w_{\mathrm{gen}}`,  
  the factors needed to scale MC to data.

- ``compute_pileup_weights.py``  
  Compute scale factors for pileup reweighting.

- ``delete_duplicate_outputs.py``  
  Identify duplicated per-event outputs from ``select_events.py``  
  and move or delete duplicates.

- ``export_hists_from_state.py``  
  Save all histograms contained in a Pepper processor state,  
  even if processing has not finished yet.

- ``generate_btag_efficiencies.py``  
  Generate a ROOT file containing efficiency histograms  
  needed for b-tagging scale factors.

- ``generate_jet_puid_efficiencies``  
  Generate a ROOT file containing efficiency histograms  
  needed for jet pile-up ID scale factors.

- ``get_bad_local_files.py``  
  Find NanoAOD files that exist in the store directory  
  but cannot be accessed (e.g. due to technical issues).

- ``hdf5_to_ttree.py``  
  Merge and convert Pepper HDF5 files to ROOT files containing TTrees.

- ``merge_hists.py``  
  Compute weighted average of two scale-factor histograms.

- ``plot_histograms.py``  
  Plot histograms in a ratio-plot style.

- ``produce_met_xy_nums.py``  
  Convert centrally provided MET-xy correction headers into JSON files.

- ``rucio_create_rules.py``  
  Create Rucio rules for all datasets specified in a Pepper config.  
  Once approved, datasets are transferred to the local site.

- ``ttbarll_dy_sf_calculate.py``  
  Calculate DY reweighting scale factors from the output  
  of ``ttbarll_dy_sf_produce.py``.

- ``ttbarll_dy_sf_produce.py``  
  Produce the numbers needed for Drell–Yan SF calculation.

- ``ttbarll_kinreco_hists_produce.py``  
  Produce histograms needed for top-quark kinematic reconstruction.

- ``ttbarll_select_events.py``  
  Run the main ``ttbarll`` analysis, producing histograms  
  and per-event data.

- ``ttbarll_trigger_sf_calculate.py``  
  Compute SFs for dileptonic ``ttbar`` triggers using cross-trigger  
  method from the output of ``ttbarll_trigger_sf_produce.py``.

- ``ttbarll_trigger_sf_produce.py``  
  Produce the numbers needed for trigger scale-factor calculation.
