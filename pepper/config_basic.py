#!/usr/bin/env python3
import os
from pepper.misc import get_run_for_year, LHCRun
import uproot
import hjson
import coffea
import coffea.jetmet_tools
from coffea import lookup_tools
from functools import partial
import correctionlib

import pepper
from pepper.scale_factors import (
    TopPtWeigter,
    PileupWeighter,
    BTagWPs,
    BTagWeighter,
    get_evaluator,
    ScaleFactors,
    MuonScaleFactor,
    JetPuIdWeighter,
    CorrLibSFs,
    LegacyJetAdapter,
    CorrectionlibAdapter,
    CorrectionlibCompoundAdapter,
    JetIdProducer,
)


class ConfigBasicPhysics(pepper.Config):
    def __init__(self, path_or_file, textparser=hjson.load, cwd="."):
        """Initialize the configuration.

        The default list of required arguments, behaviours, and  special variables
        is:

        - Required Args (inherited from Config):
            - exp_datasets
            - mc_datasets
            - mc_lumifactors
        - Behaviors (inherited from Config and extended):
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
        - Special Vars (inherited from Config):
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
        self._check_integrity_jet_corrections()
        self.behaviors.update(
            {
                "year": self._get_year,
                "crosssections": self._get_maybe_external,
                "top_pt_reweighting": self._get_top_pt_reweighting,
                "pileup_reweighting": self._get_pileup_reweighting,
                "electron_sf": self._get_scalefactors,
                "muon_sf": self._get_muonscalefactor,
                "muon_rochester": self._get_rochester_corr,
                "btag_sf": self._get_btag_sf,
                "btag_wps": self._get_btag_wps,
                "jet_ids": self._get_jet_ids,
                "jet_puid_sf": self._get_puid_sf,
                "jme_correctionlib_corrections": self._get_jme_corrections,
                "jet_correction_mc": self._get_jet_correction_legacy,
                "jet_correction_data": self._get_jet_correction_legacy_dict,
                "jet_uncertainty": partial(
                    self._get_jet_general_legacy, evaltype="junc",
                    cls=coffea.jetmet_tools.JetCorrectionUncertainty),
                "jet_resolution": partial(
                    self._get_jet_general_legacy, evaltype="jr",
                    cls=coffea.jetmet_tools.JetResolution),
                "jet_ressf": partial(
                    self._get_jet_general_legacy, evaltype="jersf",
                    cls=coffea.jetmet_tools.JetResolutionScaleFactor),
                "jet_veto_map": self._get_scalefactor,
                "MET_xy_shifts": self._get_maybe_external,
                "crosssection_uncertainty": self._get_maybe_external,
                "reco_info_file": self._get_path,
                "store": self._get_path,
                "lumimask": self._get_path,
                "rng_seed_file": self._get_path,
            }
        )

    def _check_integrity_jet_corrections(self):
        """
        Validate that either correctionlib OR legacy corrections are specified, not both.

        Raises:
            ConfigError: If both approaches are used or neither is properly configured
        """
        config = self  # readability - self is the config object
        has_correctionlib = "jme_correctionlib_corrections" in config

        legacy_keys = [
            "jet_correction_data",
            "jet_correction_mc",
            "jet_uncertainty",
            "jet_resolution",
            "jet_ressf"
        ]

        # Check if any legacy keys exist at top level and have non-empty values
        has_legacy = any(
            key in config and config[key] and (
                # Non-empty list
                (isinstance(config[key], list) and len(config[key]) > 0) or
                # Non-empty string
                (isinstance(config[key], str) and config[key].strip()) or
                # Non-empty dict
                (isinstance(config[key], dict) and len(config[key]) > 0)
            )
            for key in legacy_keys
        )

        if has_correctionlib and has_legacy:
            raise pepper.config.ConfigError(
                "Configuration error: Cannot use both 'jme_correctionlib_corrections' "
                "and legacy correction keys (jet_correction_data, jet_correction_mc, etc.) "
                "Please use only one approach."
            )

        # If using correctionlib, validate its structure
        if has_correctionlib:
            clib_config = config["jme_correctionlib_corrections"]

            if "path" not in clib_config:
                raise pepper.config.ConfigError(
                    "Configuration error: 'jme_correctionlib_corrections' must contain 'path'"
                )
            mandatory_keys = ["jet_correction_data", "jet_correction_mc", "jet_uncertainty"]
            # Validate that correction name keys are present
            missing_keys = [key for key in mandatory_keys if key not in clib_config]
            if missing_keys:
                raise pepper.config.ConfigError(
                    f"Configuration error: 'jme_correctionlib_corrections' missing keys: {missing_keys}"
                )

        return has_correctionlib

    def _get_scalefactor(self, sfpath, sysnaming={
            "central": "sf", "up": "sfup", "down": "sfdown"}):
        if not isinstance(sfpath, list) or len(sfpath) != 3:
            raise pepper.config.ConfigError(
                "Scale factors needs to be list of 3-element-lists "
                "in form of [rootfile, histname, inputs] or"
                "[jsonfile, corrname, additional_args]")
        if sfpath[0].endswith(".root"):
            with uproot.open(self._get_path(sfpath[0])) as f:
                hist = f[sfpath[1]]
            return ScaleFactors.from_hist(hist, sfpath[2])
        elif sfpath[0].endswith(".json") or sfpath[0].endswith(".json.gz"):
            return CorrLibSFs(sysnaming, self._get_path(sfpath[0]),
                              sfpath[1], sfpath[2])

    def _get_scalefactors(self, value, sysnaming={
            "central": "sf", "up": "sfup", "down": "sfdown"}):
        return [self._get_scalefactor(sfpath, sysnaming) for sfpath in value]

    def _get_muonscalefactor(self, value):
        year = self["year"]
        if ("split_muon_uncertainty" not in self
                or not self["split_muon_uncertainty"]):
            if get_run_for_year(year) == LHCRun.Run3:
                # In Run 3, central is called "nominal"
                return self._get_scalefactors(
                    value, {"central": "nominal", "up": "systup",
                            "down": "systdown"})
            else:
                # In Run 2,  central is called "sf"
                return self._get_scalefactors(
                    value, {"central": "sf", "up": "systup",
                            "down": "systdown"})

        sfs = []
        for sfpath in value:
            if not sfpath[0].endswith(".root"):
                raise NotImplementedError(
                    "Split muon uncertainties only implemented "
                    "for ROOT based SFs")
            nominal = self._get_scalefactor(sfpath)
            sfpath_stat = sfpath.copy()
            sfpath_stat[1] += "_stat"
            stat = self._get_scalefactor(sfpath_stat)
            sfpath_syst = sfpath.copy()
            sfpath_syst[1] += "_syst"
            syst = self._get_scalefactor(sfpath_syst)
            sfs.append(MuonScaleFactor(nominal, stat, syst))
        return sfs

    @staticmethod
    def _get_year(value):
        return str(value)

    @staticmethod
    def _get_top_pt_reweighting(value):
        return TopPtWeigter(**value)

    def _get_pileup_reweighting(self, value):
        if isinstance(value, str):
            # Legacy root format
            with uproot.open(self._get_path(value)) as f:
                return PileupWeighter(f)
        elif isinstance(value, list):
            sysnaming = {"central": "nominal", "up": "up", "down": "down"}
            if len(value) == 2:
                path, correction_name = value
                add_args = []
            elif len(value) == 3:
                path, correction_name, add_args = value
            else:
                raise pepper.config.ConfigError(
                    "Pileup Reweighting should be a path to a root file (legacy) "
                    "or for a correctionlib a list of [jsonfile, corrname, "
                    "additional_args], where the final element is optional.")
            path = self._get_path(path)
            return CorrLibSFs(sysnaming, path, correction_name, add_args)
        else:
            raise pepper.config.ConfigError(
                "Pileup Reweighting should be a path to a root file (legacy) "
                "or for a correctionlib a list of [jsonfile, corrname, "
                "additional_args], where the final element is optional.")

    def _get_btag_wps(self, value):
        value = self._get_maybe_external(value)
        tagger = self["btag"].split(":")[0]
        year = self["year"]
        wp_dict = {k.lower(): v for k, v in value.items()}
        if tagger in wp_dict and year in wp_dict[tagger]["wps"]:
            return BTagWPs(wp_dict, tagger, year)
        else:
            return None

    def _get_btag_sf(self, value):
        weighters = []
        tagger = self["btag"].split(":")[0]
        year = self["year"]
        method = ("fixedwp"
                  if "btag_method" not in self
                  else self["btag_method"])
        ignore_missing = (self["btag_ignoremissing"]
                          if "btag_ignoremissing" in self
                          else False)
        measure_type = ("mujets"
                        if "btag_measure_type" not in self
                        else self["btag_measure_type"])
        tagger_wps = (self["btag_wps"] if "btag_wps" in self else None)
        if method == "fixedwp" and tagger_wps is None:
            raise pepper.config.ConfigError(
                f"btag_wps not in config, or does not define wps for "
                f"{tagger} in {year}.")
        for weighter_paths in value:
            paths = [self._get_path(path) for path in weighter_paths]
            btagweighter = BTagWeighter(
                paths[0], paths[1] if len(paths) > 1 else None,
                tagger=tagger, year=year, wps=tagger_wps,
                method=method, ignore_missing=ignore_missing,
                meastype=measure_type)
            weighters.append(btagweighter)
        return weighters

    def _get_puid_sf(self, value):
        if not isinstance(value, list) or len(value) > 2:
            raise pepper.config.ConfigError(
                "jet_puid_sf should either be a one element list of the SF "
                "json (if this also contains the efficienices), or a two "
                "element list of SFs, efficiency")
        return JetPuIdWeighter(*[self._get_path(path) for path in value])

    def _get_btag_corr(self, value):
        # Dump content of hjson file to dict
        if value == "":
            return None
        elif "$DATADIR" in value:
            value = value.replace(
                "$DATADIR", self._config[self.special_vars["$DATADIR"]]
            )
            with open(value) as jf:
                data = hjson.load(jf)
            return data

    def _get_rochester_corr(self, value):
        path = self._get_path(value)
        rochester_data = lookup_tools.txt_converters.convert_rochester_file(
            path, loaduncs=False
        )
        rochester = lookup_tools.rochester_lookup.rochester_lookup(
            rochester_data
        )
        return rochester

    def _get_jme_corrections(self, value):
        """Get JME corrections.
        """
        if not isinstance(value, dict):
            msg = f"'jme_correctionlib_corrections' config should be a dict, got {type(value)}"
            raise pepper.config.ConfigError(msg)
        # Manually get the path to avoid recursion error
        # and already open the correctionlib file for speed
        try:
            correctionlib_path = self._get_path(value['path'])
            corrset = correctionlib.CorrectionSet.from_file(correctionlib_path)
        except KeyError as e:
            msg = f"Missing 'path' key in 'jme_correctionlib_corrections' config: {e}"
            raise pepper.config.ConfigError(msg) from e

        behaviours = {
            "path": self._get_path,
            "jet_uncertainty_template": lambda x: x,
            "jet_correction_mc": partial(self._get_jet_correction, correctionlib_path=corrset),
            "jet_correction_data": partial(self._get_jet_correction_dict, correctionlib_path=corrset),
            "jet_uncertainty": partial(
                self._get_jet_general,
                correctionlib_path=corrset),
            "jet_resolution": partial(
                self._get_jet_general,
                correctionlib_path=corrset),
            "jet_ressf": partial(
                self._get_jet_general,
                correctionlib_path=corrset),
        }
        updated_conf = {}
        for key, value in value.items():
            if key in behaviours:
                updated_conf[key] = behaviours[key](value)
            else:
                msg = f"Unknown JME configuration key: '{key}'."
                raise pepper.config.ConfigError(msg)
        return updated_conf

    def _get_jet_correction(self, value, correctionlib_path: os.PathLike | correctionlib.CorrectionSet):
        return CorrectionlibCompoundAdapter(correctionlib_path, value)

    def _get_jet_correction_dict(self, value, correctionlib_path: os.PathLike | correctionlib.CorrectionSet):
        if isinstance(value, dict):
            corrs = {}
            for era, val in value.items():
                corrs[era] = self._get_jet_correction(val, correctionlib_path=correctionlib_path)
            return corrs
        else:
            return self._get_jet_correction(value, correctionlib_path=correctionlib_path)

    def _get_jet_correction_legacy(self, value):
        # TODO: remove when we upgrade and only support correctionlib recipes.
        evaluators = {}
        for path in value:
            path = self._get_path(path)
            evaluators.update(get_evaluator(path, "txt", "jec"))
        fjc = coffea.jetmet_tools.FactorizedJetCorrector(**evaluators)
        return LegacyJetAdapter(fjc, 'getCorrection')

    def _get_jet_correction_legacy_dict(self, value):
        # TODO: remove when we upgrade and only support correctionlib recipes.
        if isinstance(value, dict):
            corrs = {}
            for era, val in value.items():
                corrs[era] = self._get_jet_correction_legacy(val)
            return corrs
        else:
            return self._get_jet_correction_legacy(value)

    def _get_jet_general(self, value, correctionlib_path: os.PathLike | correctionlib.CorrectionSet):
        """
        Create general jet correction objects (uncertainty, resolution, SF).
        These are SIMPLE (non-compound) corrections.
        sle
        Args:
            value: For legacy: file path
                For correctionlib: str or list[str] (first element used)
            correctionlib_path: Path to the correctionlib json file
        """
        # All non-JEC types are simple (single correction)
        if isinstance(value, list):
            template_string = self._config["jme_correctionlib_corrections"]["jet_uncertainty_template"]
            return {
                v: CorrectionlibAdapter(correctionlib_path, template_string.replace("[UNC]", v)) for v in value
                }
        else:
            return CorrectionlibAdapter(correctionlib_path, value)

    def _get_jet_general_legacy(self, value, evaltype, cls):
        # Legacy code
        # TODO: remove when we upgrade and only support correctionlib recipes.
        evaluator = get_evaluator(self._get_path(value), "txt", evaltype)
        legacy_obj = cls(**evaluator)

        # Map evaltype to method name
        method_map = {
            "junc": "getUncertainty",
            "jr": "getResolution",
            "jersf": "getScaleFactor",
        }

        method_name = method_map.get(evaltype)
        if method_name is None:
            raise pepper.config.ConfigError(f"Unknown evaltype: {evaltype}")

        return LegacyJetAdapter(legacy_obj, method_name)

    def _get_jet_ids(self, value):
        # value = [jetType, jsonfile]
        return JetIdProducer(value[0], self._get_path(value[1]))
