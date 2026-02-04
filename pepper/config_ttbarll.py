#!/usr/bin/env python3

import numpy as np
import uproot
import hjson

import pepper
from pepper.scale_factors import ScaleFactors


class ConfigTTbarLL(pepper.ConfigBasicPhysics):
    def __init__(self, path_or_file, textparser=hjson.load, cwd="."):
        """Initialize the configuration.

        - Required Args (inherited from ConfigBasicPhysics):
            - exp_datasets
            - mc_datasets
            - mc_lumifactors
        - Behaviors (inherited from ConfigBasicPhysics and extended):
            - file_blacklist: load from external file if a string is given
            - local_file_blacklist: load from external file if a string is given
            - xrootd_url_blacklist: load from external file if a string is given
            - mc_lumifactors: load from external file if a string is given
            - hists: load from external file if a string is given, parse into
              HistDefinition objects
            - year: convert to string
            - crosssections: load from external file if a string is given
            - top_pt_reweighting: parse into TopPtWeighter object
            - pileup_reweighting: parse into PileupWeighter object
            - electron_sf: parse into list of ScaleFactors objects
            - muon_sf: parse into list of MuonScaleFactor objects
            - btag_sf: parse into list of BTagWeighter objects
            - btag_wps: parse into BTagWPs object
            - jet_puid_sf: parse into JetPuIdWeighter object
            - jet_correction_mc: parse into FactorizedJetCorrector object
            - jet_correction_data: parse into dict of FactorizedJetCorrector objects
            - jet_uncertainty: parse into JetCorrectionUncertainty object
            - jet_resolution: parse into JetResolution object
            - jet_ressf: parse into JetResolutionScaleFactor object
            - MET_xy_shifts: load from external file if a string is given
            - crosssection_uncertainty: load from external file if a string is given
            - reco_info_file: resolve path
            - store: resolve path
            - lumimask: resolve path
            - rng_seed_file: resolve path
            - drellyan_sf: parse into ScaleFactors object
            - trigger_sfs: parse into dict of ScaleFactors objects
        - Special Vars (inherited from ConfigBasicPhysics):
            - $DATADIR: replaced with the value of ``datadir`` in the config
            - $CONFDIR: replaced with the directory containing the config file
            - $STOREDIR: replaced with the value of ``store`` in the config

        Parameters
        ----------
        path
            Path to the file containing the configuration
        textparser
            Callable to be used to parse the text contained in path_or_file
        cwd
            Path to use as the working directory for relative paths in the
            config. The actual working directory of the process might change
        """
        super().__init__(path_or_file, textparser, cwd)

        self.behaviors.update(
            {
                "drellyan_sf": self._get_drellyan_sf,
                "trigger_sfs": self._get_trigger_sfs,
            }
        )

    def _get_drellyan_sf(self, value):
        if isinstance(value, list):
            path, histname = value
            with uproot.open(self._get_path(path)) as f:
                hist = f[histname]
            dy_sf = ScaleFactors.from_hist(hist)
        else:
            data = self._get_maybe_external(value)
            dy_sf = ScaleFactors(
                bins=data["bins"],
                factors=np.array(data["factors"]),
                factors_up=np.array(data["factors_up"]),
                factors_down=np.array(data["factors_down"]))
        return dy_sf

    def _get_trigger_sfs(self, value):
        path, histnames = value
        ret = {}
        if len(histnames) != 3:
            raise pepper.config.ConfigError(
                "Need 3 histograms for trigger scale factors. Got "
                f"{len(histnames)}")
        with uproot.open(self._get_path(path)) as f:
            for chn, histname in zip(("is_ee", "is_em", "is_mm"), histnames):
                ret[chn] = ScaleFactors.from_hist(
                    f[histname], dimlabels=["lep1_pt", "lep2_pt"])
        return ret
