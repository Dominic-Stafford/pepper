from functools import reduce
import numpy as np
import awkward as ak
import coffea
import coffea.lumi_tools
import uproot
import logging
from dataclasses import dataclass
from typing import Optional, Tuple

import pepper
from pepper import sonnenschein
import pepper.config
from pepper.misc import get_run_for_year, LHCRun


@dataclass
class VariationArg:
    """Holds information on a systematic variation of the jet energies

    Attributes
    ----------
    name
        Name of the variation
    junc
        Tuple defining the jet energy uncertainty; first element is the
        direction ("up" or "down"), second element is the uncertainty source
    jer
        Direction into which the jet energy resolution is varied. Can also be
        "central" for no variation, but smearing still being applied
    met
        Direction into which the MET uncertainty is varied
    """
    name: Optional[str] = None
    junc: Optional[Tuple[str, str]] = None
    jer: Optional[str] = "central"
    met: Optional[str] = None


logger = logging.getLogger(__name__)


class ProcessorBasicPhysics(pepper.Processor):
    """Processor containing basic object definitions and cuts useful
       for most physics analyses."""

    config_class = pepper.ConfigBasicPhysics

    def __init__(self, config, eventdir):
        super().__init__(config, eventdir)

    def get_jetmet_variation_args(self):
        """Get a list of varations that should be done for the Jet/MET
        uncertainties"""
        ret = []
        # We split the code into two branches where code may be duplicated.
        # TODO: remove the legacy branch when we upgrade pepper to only
        # support new JME implementation.
        if "jme_correctionlib_corrections" in self.config:
            jme_config = self.config["jme_correctionlib_corrections"]

            if ("jet_resolution" not in jme_config
                    or "jet_ressf" not in jme_config):
                jer = None
            else:
                jer = "central"
            ret.append(VariationArg("UncMET_up", met="up", jer=jer))
            ret.append(VariationArg("UncMET_down", met="down", jer=jer))
            if "jet_uncertainty" in jme_config:
                junc = jme_config["jet_uncertainty"]
                levels = junc.keys()
                for source in levels:
                    if source == "jes":
                        name = "Junc_"
                    else:
                        name = f"Junc{source.replace('_', '')}_"
                    ret.append(VariationArg(
                        name + "up", junc=("up", source), jer=jer))
                    ret.append(VariationArg(
                        name + "down", junc=("down", source), jer=jer))
            if ("jet_resolution" in jme_config
                    and "jet_ressf" in jme_config):
                ret.append(VariationArg("Jer_up", jer="up"))
                ret.append(VariationArg("Jer_down", jer="down"))
        else:  # TODO: delete this branch when removing support for legacy JME
            if ("jet_resolution" not in self.config
                    or "jet_ressf" not in self.config):
                jer = None
            else:
                jer = "central"
            ret.append(VariationArg("UncMET_up", met="up", jer=jer))
            ret.append(VariationArg("UncMET_down", met="down", jer=jer))
            if "jet_uncertainty" in self.config:
                junc = self.config["jet_uncertainty"]
                if "junc_sources_to_use" in self.config:
                    levels = self.config["junc_sources_to_use"]
                else:
                    levels = junc.levels
                for source in levels:
                    if source not in junc.levels:
                        raise pepper.config.ConfigError(
                            f"Source not in jet uncertainties: {source}")
                    if source == "jes":
                        name = "Junc_"
                    else:
                        name = f"Junc{source.replace('_', '')}_"
                    ret.append(VariationArg(
                        name + "up", junc=("up", source), jer=jer))
                    ret.append(VariationArg(
                        name + "down", junc=("down", source), jer=jer))
            if "jet_resolution" in self.config and "jet_ressf" in self.config:
                ret.append(VariationArg("Jer_up", jer="up"))
                ret.append(VariationArg("Jer_down", jer="down"))
        return ret

    def get_jetmet_nominal_arg(self):
        """Get a ``VariationArg`` describing the nominal evaluation"""
        if "jme_correctionlib_corrections" in self.config:
            jme_conf = self.config["jme_correctionlib_corrections"]
        else:
            jme_conf = self.config
        if "jet_resolution" in jme_conf and "jet_ressf" in jme_conf:
            return VariationArg(None)
        else:
            return VariationArg(None, jer=None)

    def columns_to_preload(self):
        """If running with systematics, explicitly preload generator weights
        since the tracing misses those"""
        columns = super().columns_to_preload()
        if self.config["compute_systematics"]:
            return columns | \
                {"LHEScaleWeight", "LHEPdfWeight", "PSWeight"}
        else:
            return columns

    def gentop(self, data):
        """Return generator-level tops."""
        part = data["GenPart"]
        part = part[(part.genPartIdxMother != -1)
                    & part.hasFlags("isLastCopy")
                    & (abs(part.pdgId) == 6)]
        part = part[ak.argsort(part.pdgId, ascending=False)]
        return part

    def do_top_pt_reweighting(self, data):
        """Top pt reweighting according to
        https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting
        """
        pt = data["gent_lc"].pt
        weiter = self.config["top_pt_reweighting"]
        rwght = weiter(pt[:, 0], pt[:, 1])
        if weiter.sys_only:
            if self.config["compute_systematics"]:
                return np.full(len(data), True), {"Top_pt_reweighting": rwght}
            else:
                return np.full(len(data), True)
        else:
            if self.config["compute_systematics"]:
                return rwght, {"Top_pt_reweighting": 1/rwght}
            else:
                return rwght

    def do_pileup_reweighting(self, dsname, data):
        """Pileup reweighting"""
        ntrueint = data["Pileup"]["nTrueInt"]
        weighter = self.config["pileup_reweighting"]
        weight = weighter(dsname=dsname, NumTrueInteractions=ntrueint)
        if self.config["compute_systematics"]:
            # If central is zero, let up and down factors also be zero
            weight_nonzero = np.where(weight == 0, np.inf, weight)
            up = weighter(dsname=dsname, NumTrueInteractions=ntrueint, variation="up")
            down = weighter(dsname=dsname, NumTrueInteractions=ntrueint, variation="down")
            sys = {"pileup": (up / weight_nonzero, down / weight_nonzero)}
            return weight, sys
        return weight

    def add_me_uncertainties(self, dsname, selector, data):
        """Matrix-element renormalization and factorization scale"""
        # Get describtion of individual columns of this branch with
        # Events->GetBranch("LHEScaleWeight")->GetTitle() in ROOT
        if "LHEScaleWeight" not in data.fields:
            logger.warning("LHEScaleWeights missing for this sample")
            return
        lheweight = data["LHEScaleWeight"]
        if len(lheweight) == 0:
            return
        if (self.config["mc_lumifactors"] and dsname + "_LHEScaleSumw"
                in self.config["mc_lumifactors"]):
            norm = self.config["mc_lumifactors"][dsname + "_LHEScaleSumw"]
            idx = pepper.misc.get_lhe_scale_idxs(len(norm))
            selector.set_systematic(
                "MEren",
                lheweight[:, idx[0]] * abs(norm[idx[0]]),
                lheweight[:, idx[1]] * abs(norm[idx[1]]))
            selector.set_systematic(
                "MEfac",
                lheweight[:, idx[2]] * abs(norm[idx[2]]),
                lheweight[:, idx[3]] * abs(norm[idx[3]]))
        elif (not self.config["mc_lumifactors"]
              and "LHEScaleWeight" in ak.fields(data)
              and ak.num(lheweight)[0] > 0):
            idx = pepper.misc.get_lhe_scale_idxs(
                len(lheweight[0]))
            selector.set_systematic(
                "MEren", lheweight[:, idx[0]],
                lheweight[:, idx[1]], norm_post=True)
            selector.set_systematic(
                "MEfac", lheweight[:, idx[2]],
                lheweight[:, idx[3]], norm_post=True)
        else:
            logger.warning("LHEScaleWeights missing for this sample")

    def add_ps_uncertainties(self, selector, data):
        """Parton shower scale uncertainties"""
        psweight = data["PSWeight"]
        if len(psweight) == 0:
            return
        num_weights = ak.num(psweight)[0]
        if num_weights == 1:
            # NanoAOD containts one 1.0 per event in PSWeight if there are no
            # PS weights available, otherwise all counts > 1.
            return
        if num_weights == 4:
            if self.config["year"].startswith("ul"):
                # Workaround for PSWeight number changed their order in
                # NanoAODv8, meaning non-UL is unaffected
                selector.set_systematic(
                    "PSisr", psweight[:, 0], psweight[:, 2])
                selector.set_systematic(
                    "PSfsr", psweight[:, 1], psweight[:, 3])
            else:
                selector.set_systematic(
                    "PSisr", psweight[:, 2], psweight[:, 0])
                selector.set_systematic(
                    "PSfsr", psweight[:, 3], psweight[:, 1])
        else:
            logger.warning(
                "Unexpected length of the PSWeight: "
                f"{num_weights}")

    def add_pdf_uncertainties(self, dsname, selector, data):
        """Add PDF uncertainties, using the methods described here:
        https://arxiv.org/pdf/1510.03865.pdf#section.6"""
        if ("LHEPdfWeight" not in data.fields
                or ak.num(data["LHEPdfWeight"])[0] == 0):
            logger.warning("LHEPdfWeights missing for this sample")
            return
        if len(data["LHEPdfWeight"]) == 0:
            return
        if "pdf_types" not in self.config:
            logger.warning("Not processing pdfs due to missing 'pdf_types' "
                           "in config!")
            return

        pdfs = data["LHEPdfWeight"]
        pdf_doc = pdfs.__doc__
        pdf_type = None

        for LHA_ID, _type in self.config["pdf_types"].items():
            if LHA_ID in pdf_doc:
                pdf_type = _type.lower()

        split_pdf_uncs = False
        if "split_pdf_uncs" in self.config:
            split_pdf_uncs = self.config["split_pdf_uncs"]

        norm_pdf_uncs_post = False
        if ("normalize_pdf_uncs" in self.config
                and self.config["normalize_pdf_uncs"]):
            if self.config["mc_lumifactors"]:
                if dsname + "_LHEPdfSumw" not in self.config["mc_lumifactors"]:
                    raise pepper.config.ConfigError(
                        "Missing lumifactors for PDF uncertainties for dataset"
                        f" '{dsname}'. Please run compute_mc_lumifactors.py "
                        "with the '-p' option.")
                norm = self.config["mc_lumifactors"][dsname + "_LHEPdfSumw"]
                pdfs = pdfs * abs(np.array(norm)[np.newaxis, :])
            else:
                norm_pdf_uncs_post = True

        # Check if sample has alpha_s variations - currently assuming number of
        # regular variations is a multiple of 10
        if len(data) == 0:
            has_as_unc = False
        else:
            has_as_unc = (len(pdfs[0]) % 10) > 1
            # Workaround for "order of scale and pdf weights not consistent"
            # See https://twiki.cern.ch/twiki/bin/view/CMS/MCKnownIssues
            if ak.mean(pdfs[0]) < 0.6:  # approximate, see if factor 2 needed
                pdfs = ak.without_parameters(pdfs)
                pdfs = ak.concatenate([pdfs[:, 0:1], pdfs[:, 1:] * 2], axis=1)
        n_offset = -2 if has_as_unc else None

        if split_pdf_uncs:
            # Just output variations - user
            # will need to combine these for limit setting
            num_variation = len(pdfs[0]) + (n_offset or 0)
            if pdf_type == "true_hessian" or pdf_type == "hessian":
                # First element is central value - adjust all other
                # elements relative to this
                selector.set_systematic(
                    "PDF", *[pdfs[:, i] - pdfs[:, 0] + 1
                             for i in range(1, num_variation)],
                    scheme="numeric", norm_post=norm_pdf_uncs_post)
                if has_as_unc:
                    selector.set_systematic(
                        "PDFalphas",
                        pdfs[:, -1] - pdfs[:, 0] + 1,
                        pdfs[:, -2] - pdfs[:, 0] + 1,
                        norm_post=norm_pdf_uncs_post)
            elif pdf_type.startswith("mc"):
                selector.set_systematic(
                    "PDF", *[pdfs[:, i] for i in range(num_variation)],
                    scheme="numeric", norm_post=norm_pdf_uncs_post)
                if has_as_unc:
                    selector.set_systematic(
                        "PDFalphas", pdfs[:, -1], pdfs[:, -2],
                        norm_post=norm_pdf_uncs_post)
            elif pdf_type is None:
                raise pepper.config.ConfigError(
                    "PDF LHA Id not included in config. PDF docstring is: "
                    + pdf_doc)
            else:
                raise pepper.config.ConfigError(
                    f"PDF type {pdf_type} not recognised. Valid options "
                    "are 'True_Hessian', 'Hessian', 'MC' and 'MC_Gaussian'")
        else:
            if pdf_type == "true_hessian":
                # Treatment of true hessian uncertainties, for e.g. CTEQ
                # or HERA sets
                eigen_vals = ak.to_numpy(pdfs[:, 1:n_offset])
                eigen_vals = eigen_vals.reshape(
                    (eigen_vals.shape[0], eigen_vals.shape[1] // 2, 2))
                central, eigenvals = ak.broadcast_arrays(
                    pdfs[:, 0, None, None], eigen_vals)
                var_up = ak.max((eigen_vals - central), axis=2)
                var_up = ak.where(var_up > 0, var_up, 0)
                var_up = np.sqrt(ak.sum(var_up ** 2, axis=1))
                var_down = ak.max((central - eigen_vals), axis=2)
                var_down = ak.where(var_down > 0, var_down, 0)
                var_down = np.sqrt(ak.sum(var_down ** 2, axis=1))
                unc = None
            if pdf_type == "hessian":
                # Treatment of pseudo hessian uncertainties, for e.g. NNPDF
                # or pdf4LHC sets
                eigen_vals = ak.to_numpy(pdfs[:, 1:n_offset])
                variations = eigen_vals - pdfs[:, 0, None]
                unc = np.sqrt(ak.sum(variations ** 2, axis=1))
            elif pdf_type == "mc":
                # ak.sort produces an error here. Work-around:
                variations = np.sort(ak.to_numpy(pdfs[:, 1:n_offset]))
                nvar = ak.num(variations)[0]
                unc = (variations[:, int(round(0.841344746*nvar))]
                       - variations[:, int(round(0.158655254*nvar))]) / 2
            elif pdf_type == "mc_gaussian":
                mean = ak.mean(pdfs[:, 1:n_offset], axis=1)
                unc = np.sqrt((ak.sum(pdfs[:, 1:n_offset] - mean) ** 2)
                              / (ak.num(pdfs)[0] - (3 if n_offset else 1)))
            elif pdf_type is None:
                raise pepper.config.ConfigError(
                    "PDF LHA Id not included in config. PDF docstring is: "
                    + pdf_doc)
            else:
                raise pepper.config.ConfigError(
                    f"PDF type {pdf_type} not recognised. Valid options "
                    "are 'True_Hessian', 'Hessian', 'MC' and 'MC_Gaussian'")

            # Add PDF alpha_s uncertainties
            if has_as_unc:
                if ("combine_alpha_s" in self.config and
                        self.config["combine_alpha_s"]):
                    alpha_s_unc = (pdfs[:, -1] - pdfs[:, -2]) / 2
                    if unc is not None:
                        unc = np.sqrt(unc ** 2 + alpha_s_unc ** 2)
                    else:
                        var_up = np.sqrt(var_up ** 2 + alpha_s_unc ** 2)
                        var_down = np.sqrt(var_down ** 2 + alpha_s_unc ** 2)
                else:
                    selector.set_systematic(
                        "PDFalphas", pdfs[:, -1], pdfs[:, -2])
            if unc is not None:
                selector.set_systematic("PDF", 1 + unc, 1 - unc)
            else:
                selector.set_systematic("PDF", 1 + var_up, 1 - var_down)

    def add_generator_uncertainies(self, dsname, selector):
        """Add MC generator uncertainties: ME, PS and PDF"""
        data = selector.data
        self.add_me_uncertainties(dsname, selector, data)
        self.add_ps_uncertainties(selector, data)
        self.add_pdf_uncertainties(dsname, selector, data)

    def crosssection_scale(self, dsname, data):
        """Cross section uncertainties. These are values depending only on the
        data set"""
        num_events = len(data)
        lumifactors = self.config["mc_lumifactors"]
        factor = np.full(num_events, lumifactors[dsname])
        if "stitching_factors" in self.config:
            for ds, sfs in self.config["stitching_factors"].items():
                if dsname.startswith(ds):
                    factor = np.full(num_events,
                                     lumifactors[ds + "_inclusive"])
                    edges = sfs["edges"] + [np.inf]
                    for i, fac in enumerate(sfs["factors"]):
                        var = pepper.hist_defns.DataPicker(sfs["axis"])(data)
                        factor[(var >= edges[i]) & (var < edges[i + 1])] *= fac
        if (self.config["compute_systematics"]
                and not ("skip_nonshape_systematics" in self.config
                         and self.config["skip_nonshape_systematics"])
                and dsname in self.config["crosssection_uncertainty"]
                and dsname not in self.config["dataset_for_systematics"]):
            xsuncerts = self.config["crosssection_uncertainty"]
            groups = set(v[0] for v in xsuncerts.values() if v is not None)
            systematics = {}
            for group in groups:
                if xsuncerts[dsname] is None or group != xsuncerts[dsname][0]:
                    continue
                uncert = xsuncerts[dsname][1]
                systematics[group + "XS"] = (np.full(num_events, 1 + uncert),
                                             np.full(num_events, 1 - uncert))
            return factor, systematics
        else:
            return factor

    def blinding(self, is_mc, data):
        """Skip every nth event in the experimental data. One way to blind your
        analysis"""
        if not is_mc:
            return np.mod(data["event"], self.config["blinding_denom"]) == 0
        else:
            return np.full(len(data), 1/self.config["blinding_denom"])

    def good_lumimask(self, is_mc, dsname, data):
        """Keep only data events that are in the golden JSON files.
           For MC, add luminosity uncertainty"""
        if is_mc:
            weight = np.ones(len(data))
            if (self.config["compute_systematics"]
                    and not ("skip_nonshape_systematics" in self.config
                             and self.config["skip_nonshape_systematics"])):
                sys = {}
                if self.config["year"] in ("2018", "2016", "ul2018",
                                           "ul2016pre", "ul2016post"):
                    sys["lumi"] = (np.full(len(data), 1 + 0.025),
                                   np.full(len(data), 1 - 0.025))
                elif self.config["year"] in ("2017", "ul2017"):
                    sys["lumi"] = (np.full(len(data), 1 + 0.023),
                                   np.full(len(data), 1 - 0.023))
                return weight, sys
            else:
                return weight
        elif "lumimask" not in self.config:
            return np.full(len(data), True)
        else:
            run = np.array(data["run"])
            luminosity_block = np.array(data["luminosityBlock"])
            lumimask = coffea.lumi_tools.LumiMask(self.config["lumimask"])
            return lumimask(run, luminosity_block)

    def get_era(self, data, is_mc):
        """Return data-taking eras based on run number."""
        if is_mc:
            return self.config["year"] + "MC"
        else:
            if len(data) == 0:
                return "no_events"
            run = np.array(data["run"])[0]
            # Assumes all runs in file come from same era
            for era, startstop in self.config["data_eras"].items():
                if ((run >= startstop[0]) & (run <= startstop[1])):
                    return era
            raise ValueError(f"Run {run} does not correspond to any era")

    def passing_trigger(self, pos_triggers, neg_triggers, data):
        """Return mask for events that pass the trigger."""
        hlt = data["HLT"]
        available = ak.fields(hlt)
        triggered = np.full(len(data), False)
        for trigger_path in pos_triggers:
            if trigger_path not in available:
                logger.debug(f"HLT_{trigger_path} not in file")
                continue
            triggered |= np.asarray(hlt[trigger_path])
        for trigger_path in neg_triggers:
            if trigger_path not in available:
                continue
            triggered &= ~np.asarray(hlt[trigger_path])
        return triggered

    def add_l1_prefiring_weights(self, data):
        """Prefiring weights needed for 2016 and 2017 data"""
        w = data["L1PreFiringWeight"]
        nom = w["Nom"]
        if self.config["compute_systematics"]:
            sys = {"L1prefiring": (w["Up"] / nom, w["Dn"] / nom)}
            return nom, sys
        return nom

    def mpv_quality(self, data):
        """Check quality of primary vertex."""
        # Does not include check for fake. Is this even needed?
        R = np.hypot(data["PV_x"], data["PV_y"])
        return ((data["PV_chi2"] != 0)
                & (data["PV_ndof"] > 4)
                & (abs(data["PV_z"]) <= 24)
                & (R <= 2))

    def met_filters(self, is_mc, data):
        """Apply met filters."""
        year = str(self.config["year"]).lower()
        if not self.config["apply_met_filters"]:
            return np.full(data.shape, True)
        passing_filters = (
            data["Flag"]["goodVertices"]
            & data["Flag"]["globalSuperTightHalo2016Filter"]
            & data["Flag"]["EcalDeadCellTriggerPrimitiveFilter"]
            & data["Flag"]["BadPFMuonFilter"])

        if not (year in ("2016", "2017", "2018") and is_mc):
            passing_filters = (
                passing_filters & data["Flag"]["eeBadScFilter"])

        if get_run_for_year(year) == LHCRun.Run2:
            passing_filters = (
                passing_filters & data["Flag"]["HBHENoiseFilter"]
                & data["Flag"]["HBHENoiseIsoFilter"])
        elif get_run_for_year(year) == LHCRun.Run3:
            passing_filters = (
                passing_filters & data["Flag"]["BadPFMuonDzFilter"]
                & data["Flag"]["hfNoisyHitsFilter"])
        if year in ("2018", "2017"):
            passing_filters = (
                passing_filters & data["Flag"]["ecalBadCalibFilterV2"])
        if year in ("ul2018", "ul2017"):
            passing_filters = (
                passing_filters & data["Flag"]["ecalBadCalibFilter"])

        return passing_filters

    def in_transreg(self, abs_eta):
        """Check if object is in detector transition region."""
        return (1.444 < abs_eta) & (abs_eta < 1.566)

    def electron_id(self, e_id, electron):
        """Check if electrons have ID specified in the config file."""
        year = self.config["year"]
        if e_id == "skip":
            has_id = True
        elif e_id == "cut:loose":
            has_id = electron["cutBased"] >= 2
        elif e_id == "cut:medium":
            has_id = electron["cutBased"] >= 3
        elif e_id == "cut:tight":
            has_id = electron["cutBased"] >= 4
        elif e_id == "mva:noIso80":
            if get_run_for_year(year) == LHCRun.Run3:
                has_id = electron["mvaNoIso_WP80"]
            else:
                has_id = electron["mvaFall17V2noIso_WP80"]
        elif e_id == "mva:noIso90":
            if get_run_for_year(year) == LHCRun.Run3:
                has_id = electron["mvaNoIso_WP90"]
            else:
                has_id = electron["mvaFall17V2noIso_WP90"]
        elif e_id == "mva:Iso80":
            if get_run_for_year(year) == LHCRun.Run3:
                has_id = electron["mvaIso_WP80"]
            else:
                has_id = electron["mvaFall17V2Iso_WP80"]
        elif e_id == "mva:Iso90":
            if get_run_for_year(year) == LHCRun.Run3:
                has_id = electron["mvaIso_WP90"]
            else:
                has_id = electron["mvaFall17V2Iso_WP90"]
        else:
            raise ValueError("Invalid electron id string")
        return has_id

    def electron_cuts(self, electron, good_lep):
        """Apply some basic electron quality cuts.
        If ``good_lep`` is True, config values prefixed with ``good_`` for pt and ID
        are used. Otherwise the ones with prefix ``additional_`` are used."""
        if self.config["ele_cut_transreg"]:
            sc_eta_abs = abs(electron["eta"]
                             + electron["deltaEtaSC"])
            is_in_transreg = self.in_transreg(sc_eta_abs)
        else:
            is_in_transreg = np.array(False)
        if good_lep:
            e_id, pt_min = self.config[[
                "good_ele_id", "good_ele_pt_min"]]
        else:
            e_id, pt_min = self.config[[
                "additional_ele_id", "additional_ele_pt_min"]]
        eta_min, eta_max = self.config[["ele_eta_min", "ele_eta_max"]]
        return (self.electron_id(e_id, electron)
                & (~is_in_transreg)
                & (eta_min < electron["eta"])
                & (electron["eta"] < eta_max)
                & (pt_min < electron["pt"]))

    def pick_electrons(self, data):
        """Get electrons that pass basic quality cuts."""
        electrons = data["Electron"]
        return electrons[self.electron_cuts(electrons, good_lep=True)]

    def muon_id(self, m_id, muon):
        """Check if muons have ID specified in the config file."""
        if m_id == "skip":
            has_id = True
        elif m_id == "cut:loose":
            has_id = muon["looseId"]
        elif m_id == "cut:medium":
            has_id = muon["mediumId"]
        elif m_id == "cut:tight":
            has_id = muon["tightId"]
        elif m_id == "mva:loose":
            has_id = muon["mvaId"] >= 1
        elif m_id == "mva:medium":
            has_id = muon["mvaId"] >= 2
        elif m_id == "mva:tight":
            has_id = muon["mvaId"] >= 3
        else:
            raise ValueError("Invalid muon id string")
        return has_id

    def muon_iso(self, iso, muon):
        """Check if muons have isolation specified in the config file."""
        if iso == "skip":
            return True
        elif iso == "cut:very_loose":
            return muon["pfIsoId"] > 0
        elif iso == "cut:loose":
            return muon["pfIsoId"] > 1
        elif iso == "cut:medium":
            return muon["pfIsoId"] > 2
        elif iso == "cut:tight":
            return muon["pfIsoId"] > 3
        elif iso == "cut:very_tight":
            return muon["pfIsoId"] > 4
        elif iso == "cut:very_very_tight":
            return muon["pfIsoId"] > 5
        else:
            iso, iso_value = iso.split(":")
            value = float(iso_value)
            if iso == "dR<0.3_chg":
                return muon["pfRelIso03_chg"] < value
            elif iso == "dR<0.3_all":
                return muon["pfRelIso03_all"] < value
            elif iso == "dR<0.4_all":
                return muon["pfRelIso04_all"] < value
            elif iso == "miniIso":
                return muon["miniPFRelIso_all"] < value
        raise ValueError("Invalid muon iso string")

    def muon_cuts(self, muon, good_lep):
        """Apply some basic muon quality cuts
        If ``good_lep`` is True, config values prefixed with ``good_`` for pt, ID
        and iso are used. Otherwise the ones with prefix ``additional_`` are
        used."""
        if self.config["muon_cut_transreg"]:
            is_in_transreg = self.in_transreg(abs(muon["eta"]))
        else:
            is_in_transreg = np.array(False)
        if good_lep:
            m_id, pt_min, iso = self.config[[
                "good_muon_id", "good_muon_pt_min", "good_muon_iso"]]
        else:
            m_id, pt_min, iso = self.config[[
                "additional_muon_id", "additional_muon_pt_min",
                "additional_muon_iso"]]
        eta_min, eta_max = self.config[["muon_eta_min", "muon_eta_max"]]
        return (self.muon_id(m_id, muon)
                & self.muon_iso(iso, muon)
                & (~is_in_transreg)
                & (eta_min < muon["eta"])
                & (muon["eta"] < eta_max)
                & (pt_min < muon["pt"]))

    def pick_muons(self, data):
        """Get muons that pass basic quality cuts."""
        muons = data["Muon"]
        return muons[self.muon_cuts(muons, good_lep=True)]

    def apply_rochester_corr(self, rng, is_mc, data):
        """Apply Rochester corrections for muons."""
        muons = data["Muon"]
        if not is_mc:
            dtSF = self.config["muon_rochester"].kScaleDT(
                muons["charge"], muons["pt"], muons["eta"], muons["phi"]
            )
            muons["pt"] = muons["pt"] * dtSF
        else:
            # if reco pt has corresponding gen pt
            mcSF1 = self.config["muon_rochester"].kSpreadMC(
                muons["charge"],
                muons["pt"],
                muons["eta"],
                muons["phi"],
                muons.matched_gen.pt,
            )
            # if reco pt has no corresponding gen pt
            counts = ak.num(muons["pt"])
            mc_rand = rng.uniform(size=ak.sum(counts))
            mc_rand = ak.unflatten(mc_rand, counts)
            mcSF2 = self.config["muon_rochester"].kSmearMC(
                muons["charge"],
                muons["pt"],
                muons["eta"],
                muons["phi"],
                muons["nTrackerLayers"],
                mc_rand,
            )
            # Combine the two scale factors and scale the pt
            mcSF = ak.where(
                ak.is_none(muons.matched_gen.pt, axis=1), mcSF2, mcSF1)
            # Remove masking from layout, none of the SF are masked here
            mcSF = ak.fill_none(mcSF, 0)
            muons["pt"] = muons["pt"] * mcSF
        return muons

    def build_lepton_column(self, data):
        """Build a lepton column containing electrons and muons."""
        electron = data["Electron"]
        muon = data["Muon"]
        columns = ["pt", "eta", "phi", "mass", "pdgId"]
        lepton = {}
        for column in columns:
            lepton[column] = ak.concatenate([electron[column], muon[column]],
                                            axis=1)
        lepton = ak.zip(lepton, with_name="PtEtaPhiMLorentzVector",
                        behavior=data.behavior)

        # Sort leptons by pt
        # Also workaround for awkward bug using ak.values_astype
        # https://github.com/scikit-hep/awkward-1.0/issues/1288
        lepton = lepton[
            ak.values_astype(ak.argsort(lepton["pt"], ascending=False), int)]
        return lepton

    def compute_lepton_sf(self, data):
        """Compute identification and isolation scale factors for
           leptons (electrons and muons)."""
        eles = data["Electron"]
        muons = data["Muon"]

        weight = np.ones(len(data))
        systematics = {}
        # Electron identification efficiency
        for i, sffunc in enumerate(self.config["electron_sf"]):
            sceta = eles.eta + eles.deltaEtaSC
            params = {}
            for dimlabel in sffunc.dimlabels:
                if dimlabel == "abseta":
                    params["abseta"] = abs(sceta)
                elif dimlabel == "eta":
                    params["eta"] = sceta
                else:
                    params[dimlabel] = getattr(eles, dimlabel)
            central = ak.prod(sffunc(**params), axis=1)
            key = "electronsf{}".format(i)
            if self.config["compute_systematics"]:
                up = ak.prod(sffunc(**params, variation="up"), axis=1)
                down = ak.prod(sffunc(**params, variation="down"), axis=1)
                systematics[key] = (up / central, down / central)
            weight = weight * central
        # Muon identification and isolation efficiency
        for i, sffunc in enumerate(self.config["muon_sf"]):
            params = {}
            for dimlabel in sffunc.dimlabels:
                if dimlabel == "abseta":
                    params["abseta"] = abs(muons.eta)
                else:
                    params[dimlabel] = getattr(muons, dimlabel)
            central = ak.prod(sffunc(**params), axis=1)
            key = f"muonsf{i}"
            if self.config["compute_systematics"]:
                if ("split_muon_uncertainty" not in self.config
                        or not self.config["split_muon_uncertainty"]):
                    unctypes = ("",)
                else:
                    unctypes = ("stat ", "syst ")
                for unctype in unctypes:
                    up = ak.prod(sffunc(
                        **params, variation=f"{unctype}up"), axis=1)
                    down = ak.prod(sffunc(
                        **params, variation=f"{unctype}down"), axis=1)
                    systematics[key + unctype.replace(" ", "")] = (
                        up / central, down / central)
            weight = weight * central
        return weight, systematics

    def lepton_pair(self, is_mc, data):
        """Select events that contain at least two leptons."""
        accept = np.asarray(ak.num(data["Lepton"]) >= 2)
        if is_mc:
            weight, systematics = self.compute_lepton_sf(data[accept])
            accept = accept.astype(float)
            accept[accept.astype(bool)] *= np.asarray(weight)
            return accept, systematics
        else:
            return accept

    def opposite_sign_lepton_pair(self, data):
        """Select events that contain two opposite-sign leptons."""
        return (np.sign(data["Lepton"][:, 0].pdgId)
                != np.sign(data["Lepton"][:, 1].pdgId))

    def same_flavor_lepton_pair(self, data):
        """Select events that contain two same-flavor leptons."""
        return (abs(data["Lepton"][:, 0].pdgId)
                == abs(data["Lepton"][:, 1].pdgId))

    def mass_lepton_pair(self, data):
        """Return invariant mass of lepton pair."""
        return (data["Lepton"][:, 0] + data["Lepton"][:, 1]).mass

    def compute_jec_factor(self, is_mc, era, name_jet, data, pt=None, eta=None,
                           phi=None, area=None, rho=None, raw_factor=None):
        """Return jet energy correction factor."""
        if name_jet != "Jet" and name_jet != "FatJet":
            raise ValueError("{} is not allowed for jec. Viable options: Jet, FatJet"
                             .format(name_jet))
        if pt is None:
            pt = data[name_jet].pt
        if eta is None:
            eta = data[name_jet].eta
        if phi is None:
            phi = data[name_jet].phi
        if area is None:
            area = data[name_jet].area
        if rho is None:
            if "Rho" in data.fields:
                rho = data["Rho"]["fixedGridRhoFastjetAll"]
            else:
                rho = data["fixedGridRhoFastjetAll"]
        if raw_factor is None:
            raw_factor = 1 - data[name_jet]["rawFactor"]
        if is_mc:
            if name_jet == "Jet":
                if self.config.get('jme_correctionlib_corrections'):
                    jec = self.config['jme_correctionlib_corrections']["jet_correction_mc"]
                else:
                    jec = self.config["jet_correction_mc"]
            elif name_jet == "FatJet":  # Not `else` for readability
                if self.config.get('fatjet_jme_correctionlib_corrections'):
                    jec = self.config['fatjet_jme_correctionlib_corrections']["jet_correction_mc"]
                else:
                    raise ValueError("FatJet JEC requires 'fatjet_jme_correctionlib_corrections'!")
        else:
            if name_jet == "Jet":
                if self.config.get('jme_correctionlib_corrections'):
                    jec = self.config['jme_correctionlib_corrections']["jet_correction_data"]
                else:
                    jec = self.config["jet_correction_data"]
            elif name_jet == "FatJet":  # Not `else` for readability
                if self.config.get('fatjet_jme_correctionlib_corrections'):
                    jec = self.config['fatjet_jme_correctionlib_corrections']["jet_correction_data"]
                else:
                    raise ValueError("FatJet JEC requires 'fatjet_jme_correctionlib_corrections'!")
            if isinstance(jec, dict):
                jec = jec[era]

        raw_pt = pt * raw_factor
        if (self.config.get("jme_correctionlib_corrections")
                or self.config.get('fatjet_jme_correctionlib_corrections')):
            l1l2l3 = jec(
                JetPt=raw_pt, JetEta=eta, JetPhi=phi, JetA=area, Rho=rho, run=data.run)
        elif self.config["year"] == "2023post":
            l1l2l3 = jec(
                JetPt=raw_pt, JetEta=eta, JetPhi=phi, JetA=area, Rho=rho)
        else:
            l1l2l3 = jec(
                JetPt=raw_pt, JetEta=eta, JetA=area, Rho=rho)
        return raw_factor * l1l2l3

    def get_junc_factor_mask(self, data, source, pt, eta, flavor):
        """
        Some Jet enery uncertainties are for specific jet flavors only. This
        method decides which jets to apply the uncertainties to.

        Returns None to indicate to apply the uncertainty to all jets.
        Otherwise a bool ak.Array with either True or False for every jet
        """
        return None  # To be modified in subclasses

    def compute_junc_factor(self, data, name_jet, variation, source="jes", pt=None,
                            eta=None, flavor=None):
        """Return jet energy correction uncertainty factor."""
        if name_jet != "Jet" and name_jet != "FatJet":
            raise ValueError("{} is not allowed for junc. Viable options: Jet, FatJet"
                             .format(name_jet))
        if variation not in ("up", "down"):
            raise ValueError("variation must be either 'up' or 'down'")
        if pt is None:
            pt = data[name_jet].pt
        if eta is None:
            eta = data[name_jet].eta
        if (name_jet == "Jet" and flavor is None):  # FatJet doesn't support parton flavor
            flavor = data[name_jet].partonFlavour
        counts = ak.num(pt)
        if ak.sum(counts) == 0:
            return ak.unflatten([], counts)
        if name_jet == "Jet":
            if "jme_correctionlib_corrections" in self.config:
                junc = self.config['jme_correctionlib_corrections']["jet_uncertainty"][source](
                        JetPt=pt, JetEta=eta)
                # Symmetric according to this:
                # https://cms-jerc.web.cern.ch/JECUncertaintySources/#description
                junc = 1+junc if variation == "up" else 1-junc
            else:
                junc = dict(self.config["jet_uncertainty"](
                    JetPt=pt, JetEta=eta))[source]
                junc = junc[:, :, 0 if variation == "up" else 1]
            mask = self.get_junc_factor_mask(data, source, pt, eta, flavor)
            if mask is not None:
                junc = ak.where(mask, junc, ak.ones_like(junc))
        elif name_jet == "FatJet":  # Not `else` for readability
            if "fatjet_jme_correctionlib_corrections" in self.config:
                junc = self.config['fatjet_jme_correctionlib_corrections']["jet_uncertainty"][source](
                        JetPt=pt, JetEta=eta)
                junc = 1+junc if variation == "up" else 1-junc
            else:
                raise ValueError("FatJet JUNC requires 'fatjet_jme_correctionlib_corrections'!")
        return junc

    def find_matched_genjet(self, jer, jets, r=0.4):
        """Find a matched GenJet for the purpose of JER smearing,
        passing the requirements for matching according to JME.
        See https://cms-jerc.web.cern.ch/JER """
        genjet = jets.matched_gen
        deltar = jets.delta_r(genjet)
        rel_dpt = abs(jets.pt - genjet.pt)/jets.pt

        is_matched = (deltar < 0.5*r) & (rel_dpt < 3 * jer)
        genjets_matched = ak.mask(genjet, is_matched)
        return genjets_matched

    def compute_jer_factor(self, data, rng, name_jet, variation="central", pt=None,
                           eta=None, hybrid=True):
        """Return jet energy resolution factor."""
        # Coffea offers a class named JetTransformer for this. Unfortunately
        # it is more inconvinient and bugged than useful.
        if name_jet != "Jet" and name_jet != "FatJet":
            raise ValueError("{} is not allowed for jer. Viable options: Jet, FatJet"
                             .format(name_jet))
        if pt is None:
            pt = data[name_jet].pt
        if eta is None:
            eta = data[name_jet].eta
        counts = ak.num(pt)
        if ak.sum(counts) == 0:
            return ak.unflatten([], counts)
        if "Rho" in data.fields:
            rho = data["Rho"]["fixedGridRhoFastjetAll"]
        else:
            rho = data["fixedGridRhoFastjetAll"]
        if name_jet == "Jet":
            if self.config.get('jme_correctionlib_corrections'):
                jer = self.config["jme_correctionlib_corrections"]["jet_resolution"](
                    JetPt=pt, JetEta=eta, Rho=rho)
                # this handles the JER unc in the JEC format:
                # https://cms-jerc.web.cern.ch/JER/#smearing-procedures
                if "jet_ressf_uncertainty" in self.config["jme_correctionlib_corrections"]:
                    jersf = self.config["jme_correctionlib_corrections"]["jet_ressf"](
                        JetPt=pt, JetEta=eta)
                    jersf_unc = self.config["jme_correctionlib_corrections"]["jet_ressf_uncertainty"](
                        JetPt=pt, JetEta=eta)
                    if variation == "central":
                        jersf = jersf
                    elif variation == "up":
                        jersf = jersf*(1+jersf_unc)
                    elif variation == "down":
                        jersf = jersf*(1-jersf_unc)
                    else:
                        raise ValueError("variation must be one of 'central', 'up' or 'down'")
                else:
                    jersf = self.config["jme_correctionlib_corrections"]["jet_ressf"](
                        JetPt=pt, JetEta=eta, Rho=rho, systematic=variation)
            else:
                jer = self.config["jet_resolution"](
                    JetPt=pt, JetEta=eta, Rho=rho)
                jersf = self.config["jet_ressf"](
                    JetPt=pt, JetEta=eta, Rho=rho)
                if variation == "central":
                    jersf = jersf[:, :, 0]
                elif variation == "up":
                    jersf = jersf[:, :, 1]
                elif variation == "down":
                    jersf = jersf[:, :, 2]
                else:
                    raise ValueError("variation must be one of 'central', 'up' or 'down'")
        elif name_jet == "FatJet":  # Not `else` for readability
            if self.config.get('fatjet_jme_correctionlib_corrections'):
                jer = self.config["fatjet_jme_correctionlib_corrections"]["jet_resolution"](
                        JetPt=pt, JetEta=eta, Rho=rho)
                jersf = self.config["fatjet_jme_correctionlib_corrections"]["jet_ressf"](
                        JetPt=pt, JetEta=eta, Rho=rho, systematic=variation)
            else:
                raise ValueError("FatJet JER requires 'fatjet_jme_correctionlib_corrections'!")
        jersmear = jer * rng.normal(size=len(jer))
        factor_stoch = 1 + np.sqrt(np.maximum(jersf**2 - 1, 0)) * jersmear
        if hybrid:
            # Hybrid method: Apply scaling relative to genpt if possible
            if name_jet == "Jet":
                matched_genjets = self.find_matched_genjet(jer, data[name_jet], 0.4)
            elif name_jet == "FatJet":  # Not `else` for readability
                matched_genjets = self.find_matched_genjet(jer, data[name_jet], 0.8)
            genpt = matched_genjets.pt
            factor_scale = 1 + (jersf - 1) * (pt - genpt) / pt
            factor = ak.where(
                ak.is_none(genpt, axis=1), factor_stoch, factor_scale)
            factor = ak.drop_none(factor)
        else:
            factor = factor_stoch
        return factor

    def compute_jet_factors(self, is_mc, era, jec, junc, jer, rng, data):
        """Return total jet factor."""
        factor = ak.ones_like(data["Jet"].pt)
        if jec:
            jecfac = self.compute_jec_factor(is_mc, era, "Jet", data)
            factor = factor * jecfac
        if is_mc and junc is not None:
            juncfac = self.compute_junc_factor(data, "Jet", *junc)
            factor = factor * juncfac
        if is_mc and jer is not None:
            jerfac = self.compute_jer_factor(data, rng, "Jet", jer)
            factor = factor * jerfac
        ret = {}
        if jec or (is_mc and (junc is not None or jer is not None)):
            ret["jetfac"] = factor
        if is_mc and jer is not None:
            ret["jerfac"] = jerfac
        return ret

    def compute_fatjet_factors(self, is_mc, era, jec, junc, jer, rng, data):
        """Return total fatjet factor."""
        fatjet_factor = ak.ones_like(data["FatJet"].pt)
        if jec:
            jecfac = self.compute_jec_factor(is_mc, era, "FatJet", data)
            fatjet_factor = fatjet_factor * jecfac
        if is_mc and junc is not None:
            juncfac = self.compute_junc_factor(data, "FatJet", *junc)
            fatjet_factor = fatjet_factor * juncfac
        if is_mc and jer is not None:
            fatjet_jerfac = self.compute_jer_factor(data, rng, "FatJet", jer)
            fatjet_factor = fatjet_factor * fatjet_jerfac
        ret = {}
        if jec or (is_mc and (junc is not None or jer is not None)):
            ret["fatjetfac"] = fatjet_factor
        if is_mc and jer is not None:
            ret["fatjerfac"] = fatjet_jerfac
        return ret

    def _evaluate_jet_ids(self, data, collection, config_key, prefix):
        """Evaluate the jet IDs, which are external correctionlib files
        starting from nano v15. They are set as separate columns since
        they are needed for all jets (not just the ones passing the
        analysis requirements) for the jet veto maps."""
        jets = data[collection]
        if "jetId" in jets.fields:  # Legacy jet ID
            return {
                f"{prefix}Loose": jets.isLoose,
                f"{prefix}Tight": jets.isTight,
                f"{prefix}TightLeptonVeto": jets.isTightLeptonVeto,
            }
        else:
            jet_id_evaluator = self.config[config_key]
            tight_id, tight_lep_veto_id = jet_id_evaluator.evaluate(jets)
            return {
                f"{prefix}Tight": tight_id,
                f"{prefix}TightLeptonVeto": tight_lep_veto_id,
            }

    def evaluate_jet_ids(self, data):
        """Evaluate jet ID for R=0.4 jets"""
        return self._evaluate_jet_ids(
            data,
            collection="Jet",
            config_key="jet_ids",
            prefix="JetId",
        )

    def evaluate_fatjet_ids(self, data):
        """Evaluate jet ID for R=0.8 jets"""
        return self._evaluate_jet_ids(
            data,
            collection="FatJet",
            config_key="fatjet_ids",
            prefix="FatJetId",
        )

    def apply_jet_veto_map(self, data):
        """Apply jet veto map according to the recommendations here:
        https://cms-jerc.web.cern.ch/Recommendations/#jet-veto-maps.
        First, a loose nominal selection is made. Then, veto events if any
        jet fulfills the loose selection criteria and lies in the jet veto
        region."""
        jets = data["Jet"]
        jetid = data["JetIdTightLeptonVeto"]
        nominal_selection_mask = (
            (jets.pt > 15)
            & jetid
            & ((jets.chEmEF + jets.neEmEF) < 0.9)
        )
        veto_map = self.config["jet_veto_map"]
        veto_jets = veto_map(eta=jets.eta[nominal_selection_mask],
                             phi=jets.phi[nominal_selection_mask])
        veto_event = ak.any(veto_jets, axis=1)
        return ~veto_event

    def _good_jet_like(self, data, collection, jet_key, pt_factor_name=None):
        """Apply some basic quality cuts to a jet-like collection."""
        jets = data[collection]
        leptons = data["Lepton"]
        j_id, lep_dist, eta_min, eta_max, pt_min = self.config[[
            f"good_{jet_key}_id",
            f"good_{jet_key}_lepton_distance",
            f"good_{jet_key}_eta_min",
            f"good_{jet_key}_eta_max",
            f"good_{jet_key}_pt_min",
        ]]

        # Jet IDs should be set in the selector using evaluate_jet_ids or evaluate_fatjet_ids
        if j_id == "skip":
            has_id = True
        elif j_id == "cut:loose":
            if f"{collection}IdLoose" not in data.fields:
                raise ValueError(f"Loose {collection} ID is not supported for this data era.")
            has_id = data[f"{collection}IdLoose"]
            # Always False in 2017 and 2018
        elif j_id == "cut:tight":
            has_id = data[f"{collection}IdTight"]
        elif j_id == "cut:tightlepveto":
            has_id = data[f"{collection}IdTightLeptonVeto"]
        else:
            raise pepper.config.ConfigError(
                "Invalid good_jet_id: {}".format(j_id))

        j_pt = jets.pt
        if pt_factor_name is not None and pt_factor_name in ak.fields(data):
            j_pt = j_pt * data[pt_factor_name]
        has_lepton_close = ak.any(
            jets.metric_table(leptons) < lep_dist, axis=2)

        return (has_id
                & (~has_lepton_close)
                & (eta_min < jets.eta)
                & (jets.eta < eta_max)
                & (pt_min < j_pt))

    def good_jet(self, data):
        """Apply jet quality cuts to R=0.4 jets."""
        return self._good_jet_like(
            data,
            collection="Jet",
            jet_key="jet",
            pt_factor_name="jetfac"
        )

    def good_fatjet(self, data):
        """Apply jet quality cuts to R=0.8 jets."""
        return self._good_jet_like(
            data,
            collection="FatJet",
            jet_key="fatjet",
            pt_factor_name="fatjetfac"
        )

    def has_puid(self, jets):
        """Whether jets satisfy the configured pileup ID"""
        j_puId = self.config["good_jet_puId"]
        is2016 = ("ul2016" in self.config["year"])
        if j_puId == "skip":
            has_puId = True
        elif ((not is2016 and j_puId == "cut:loose")
              or (is2016 and j_puId == "cut:tight")):
            has_puId = ak.values_astype(jets["puId"] & 0b100, bool)
        elif j_puId == "cut:medium":
            has_puId = ak.values_astype(jets["puId"] & 0b10, bool)
        elif ((not is2016 and j_puId == "cut:tight")
              or (is2016 and j_puId == "cut:loose")):
            has_puId = ak.values_astype(jets["puId"] & 0b1, bool)
        else:
            raise pepper.config.ConfigError(
                "Invalid good_jet_id: {}".format(j_puId))
        # Only apply PUID if pT < 50 GeV
        has_puId = has_puId | (jets.pt >= 50)
        return has_puId

    def build_jet_column(self, is_mc, data):
        """Build a column of jets passing the jet quality cuts,
           including a 'btag' key (containing the value of the
           chosen btag discriminator) and a 'btagged' key
           (to select jets that are tagged as b-jets)."""
        is_good_jet = self.good_jet(data)
        jets = data["Jet"][is_good_jet]
        if "jetfac" in ak.fields(data):
            jets["pt"] = jets["pt"] * data["jetfac"][is_good_jet]
            jets["mass"] = jets["mass"] * data["jetfac"][is_good_jet]
            jets = jets[ak.argsort(jets["pt"], ascending=False)]

        if not self.config["btag"] == "skip":
            # Evaluate b-tagging
            tagger, wp = self.config["btag"].split(":")
            if tagger == "deepcsv":
                jets["btag"] = jets["btagDeepB"]
            elif tagger == "deepjet":
                jets["btag"] = jets["btagDeepFlavB"]
            elif tagger == "robustparticletransformer":
                jets["btag"] = jets["btagRobustParTAK4B"]
            elif tagger == "particlenet":
                jets["btag"] = jets["btagPNetB"]
            elif tagger == "upart":
                jets["btag"] = jets["btagUParTAK4B"]
            else:
                raise pepper.config.ConfigError(
                    "Invalid tagger name: {}".format(tagger))
            year = self.config["year"]
            if "btag_wps" in self.config and self.config["btag_wps"]:
                _, _, wp_val = self.config["btag_wps"](wp)
            else:
                raise pepper.config.ConfigError(
                    f"btag_wps not in config, or does not define wps for "
                    f"{tagger} in {year}.")
            jets["btagged"] = jets["btag"] > wp_val
        jets["pass_pu_id"] = self.has_puid(jets)
        if is_mc:
            # A jet is considered to be a pileup jet if there is no gen jet
            # within Delta R < 0.4
            jets["has_gen_jet"] = ak.fill_none(
                jets.delta_r(jets.matched_gen) < 0.4, False)
        return jets

    def jets_with_puid(self, data):
        """Get all jets satisfying the configured pileup ID"""
        jets = data["Jet"]
        return jets[jets.pass_pu_id]

    def contain_qqb(self, fatjet, part):
        """Match fat jets to qqb from top decay"""
        # Parton is from the hard-process before radiation
        part = part[part.hasFlags("isFirstCopy", "fromHardProcess")]
        part = part[part.genPartIdxMother >= 0]  # Parton is not from the initial state
        mother = part.distinctParent
        grandmom = mother.distinctParent
        partId = abs(part.pdgId)
        momId = abs(mother.pdgId)
        # For ttW, the W might be offshell and the mother becomes an initial parton
        grandId = ak.fill_none(abs(grandmom.pdgId), 0)
        b_from_top = part[(partId == 5) & (momId == 6)]
        q_frW_frT = part[(partId < 6) & (momId == 24) & (grandId == 6)]

        dR_B = ak.fill_none(fatjet.metric_table(b_from_top), 10.)
        B_inside = ak.any(dR_B < 0.8, axis=2)

        dR_Qt = ak.fill_none(fatjet.metric_table(q_frW_frT), 10.)
        num_Qt_inside = ak.num(dR_Qt[dR_Qt < 0.8], axis=2)
        Qt_inside = (num_Qt_inside >= 2)

        return (B_inside & Qt_inside)

    def contain_qq(self, fatjet, part):
        """Match fat jets to qq from W decay or top jets w/o b inside"""
        # Parton is from the hard-process before radiation
        part = part[part.hasFlags("isFirstCopy", "fromHardProcess")]
        part = part[part.genPartIdxMother >= 0]  # Parton is not from the initial state
        mother = part.distinctParent
        grandmom = mother.distinctParent
        partId = abs(part.pdgId)
        momId = abs(mother.pdgId)
        # For ttW, the W might be offshell and the mother becomes an initial parton
        grandId = ak.fill_none(abs(grandmom.pdgId), 0)
        b_from_top = part[(partId == 5) & (momId == 6)]
        q_frW_frT = part[(partId < 6) & (momId == 24) & (grandId == 6)]
        q_frW = part[(partId < 6) & (momId == 24)]

        dR_B = ak.fill_none(fatjet.metric_table(b_from_top), 10.)
        B_outside = ak.all(dR_B >= 0.8, axis=2)

        dR_Qt = ak.fill_none(fatjet.metric_table(q_frW_frT), 10.)
        num_Qt_inside = ak.num(dR_Qt[dR_Qt < 0.8], axis=2)
        Qt_inside = (num_Qt_inside >= 2)

        dR_Qw = ak.fill_none(fatjet.metric_table(q_frW), 10.)
        num_Qw_inside = ak.num(dR_Qw[dR_Qw < 0.8], axis=2)
        Qw_inside = (num_Qw_inside >= 2)

        noHadTop = (ak.num(q_frW_frT, axis=1) == 0)

        return ((B_outside & Qt_inside) | (noHadTop & Qw_inside))

    def build_fatjet_column(self, is_mc, data):
        """Build a column of AK8 jets"""
        is_good_jet = self.good_fatjet(data)
        jets = data["FatJet"][is_good_jet]
        if "fatjetfac" in ak.fields(data):
            jets["pt"] = jets["pt"] * data["fatjetfac"][is_good_jet]
            # FatJet doesn't apply JEC to jet mass
            jets = jets[ak.argsort(jets["pt"], ascending=False)]

        if is_mc:
            has_gen_close = ak.fill_none(
                    jets.delta_r(jets.matched_gen) < 0.8, False)
            jets["has_gen_fatjet"] = has_gen_close
        return jets

    def build_lowptjet_column(self, is_mc, era, junc, jer, rng, data):
        """Build a column of low-pt jets, needed to propagate jet
           corrections to low-pt jets and consequently build the
           MET column."""
        jets = data["CorrT1METJet"]
        # For MET we care about jets close to 15 GeV. JEC is derived with
        # pt > 10 GeV and |eta| < 5.2, thus cut there
        jets = jets[(jets.rawPt > 10) & (abs(jets.eta) < 5.2)]
        l1l2l3 = self.compute_jec_factor(
            is_mc, era, "Jet", data, jets.rawPt, jets.eta, jets.phi,
            jets.area, raw_factor=ak.ones_like(jets.rawPt))
        jets["pt"] = l1l2l3 * jets.rawPt
        jets["pt_nomuon"] = jets["pt"] * (1 - jets["muonSubtrFactor"])

        # Actually if reapply_met is True one would need to take into account
        # the difference between JEC applied in NanoAOD and the one that is
        # applied in the Processor for MET type-1 corrections. However, as
        # NanoAOD doesn't provide a raw factor for low pt jets, assume the
        # difference is negligible.
        if junc is not None:
            jets["juncfac"] = self.compute_junc_factor(
                data, "Jet", *junc, pt=jets["pt"], eta=jets["eta"],
                flavor=np.zeros_like(jets["pt"]))
        else:
            jets["juncfac"] = ak.ones_like(jets["pt"])
        if jer is not None:
            jets["jerfac"] = self.compute_jer_factor(
                data, rng, "Jet", jer, jets["pt"], jets["eta"], False)
        else:
            jets["jerfac"] = ak.ones_like(jets["pt"])

        jets["emef"] = jets["mass"] = ak.zeros_like(jets["pt"])
        jets.behavior = data["Jet"].behavior
        jets = ak.with_name(jets, "PtEtaPhiMLorentzVector")

        return jets

    def build_met_column(self, is_mc, junc, jer, rng, era, data,
                         variation="central"):
        """Build a column for missing transverse energy.
        If varation is 'up' or 'down', the unclustered MET varation (up or
        down) will be applied. If it is None or 'central', nominal MET is
        used."""
        nano_met_name = None
        if get_run_for_year(self.config["year"]) == LHCRun.Run3:
            default_nano_met_name = "PuppiMET"
        else:
            default_nano_met_name = "MET"
        nano_met_name = self.config.get("jet_type", default_nano_met_name)
        met = data[nano_met_name]
        metx = met.pt * np.cos(met.phi)
        mety = met.pt * np.sin(met.phi)
        if ("MET_xy_shifts" in self.config and era != "no_events"):
            metshifts = self.config["MET_xy_shifts"]
            metx = metx - (metshifts["METxcorr"][era][0]
                           * data["PV"]["npvs"]
                           + metshifts["METxcorr"][era][1])
            mety = mety - (metshifts["METycorr"][era][0]
                           * data["PV"]["npvs"]
                           + metshifts["METycorr"][era][1])
        if variation == "up":
            if nano_met_name == "PuppiMET":
                # Sometimes nans are in the variations as documented here:
                # https://cms-talk.web.cern.ch/t/nan-values-in-puppimet-variations/127652/2
                # which according to JetMET POG may be related to this issue
                # https://github.com/cms-sw/cmssw/issues/39110. Currently,
                # we just fall back to nominal MET in such cases but emit a
                # warning. TODO: Remove this workaround once the issue is fixed.
                if (
                    np.any(np.isnan(met.ptUnclusteredUp))
                    or np.any(np.isnan(met.phiUnclusteredUp))
                ):
                    import warnings
                    warnings.warn(
                        "NaN values found in PuppiMET unclustered energy "
                        "up variation. Falling back to nominal MET for "
                        "these events.")
                    met["pt"] = np.where(np.isnan(met.ptUnclusteredUp),
                                         met.pt,
                                         met.ptUnclusteredUp
                                         )
                    met["phi"] = np.where(np.isnan(met.phiUnclusteredUp),
                                          met.phi,
                                          met.phiUnclusteredUp
                                          )
                else:
                    met["pt"] = met.ptUnclusteredUp
                    met["phi"] = met.phiUnclusteredUp
                return met
            metx = metx + met.MetUnclustEnUpDeltaX
            mety = mety + met.MetUnclustEnUpDeltaY
        elif variation == "down":
            if nano_met_name == "PuppiMET":
                # See comment in 'up' variation
                if (
                    np.any(np.isnan(met.ptUnclusteredDown))
                    or np.any(np.isnan(met.phiUnclusteredDown))
                ):
                    import warnings
                    warnings.warn(
                        "NaN values found in PuppiMET unclustered energy "
                        "down variation. Falling back to nominal MET for "
                        "these events.")
                    met["pt"] = np.where(np.isnan(met.ptUnclusteredDown),
                                         met.pt,
                                         met.ptUnclusteredDown
                                         )
                    met["phi"] = np.where(np.isnan(met.phiUnclusteredDown),
                                          met.phi,
                                          met.phiUnclusteredDown
                                          )
                else:
                    met["pt"] = met.ptUnclusteredDown
                    met["phi"] = met.phiUnclusteredDown
                # TODO: reconsider this method. If PuppiMET is renamed to MET
                # in the future, this will break.
                # TODO: Implement JET corrections for PuppiMET
                met["pt"] = np.where(np.isnan(met.ptUnclusteredDown), met.pt, met.ptUnclusteredDown)
                met["phi"] = np.where(np.isnan(met.phiUnclusteredDown), met.phi, met.phiUnclusteredDown)
                return met
            metx = metx - met.MetUnclustEnUpDeltaX
            mety = mety - met.MetUnclustEnUpDeltaY
        elif variation != "central" and variation is not None:
            raise ValueError(
                "variation must be one of None, 'central', 'up' or 'down'")
        if "jetfac" in ak.fields(data):
            factors = data["jetfac"]
            if "smear_met" in self.config:
                smear_met = self.config["smear_met"]
            else:
                smear_met = False
            if "jerfac" in ak.fields(data) and not smear_met:
                factors = factors / data["jerfac"]
            # Do MET type-1 and smearing corrections
            jets = data["OrigJet"]
            jets = ak.zip({
                "pt": jets.pt,
                "pt_nomuon": jets.pt * (1 - jets.muonSubtrFactor),
                "factor": factors - 1,
                "eta": jets.eta,
                "phi": jets.phi,
                "mass": jets.mass,
                "emef": jets.neEmEF + jets.chEmEF
            }, with_name="PtEtaPhiMLorentzVector", behavior=jets.behavior)
            lowptjets = self.build_lowptjet_column(
                is_mc, era, junc, jer, rng, data)
            # Cut according to MissingETRun2Corrections Twiki
            jets = jets[(jets["pt_nomuon"] > 15) & (jets["emef"] < 0.9)]
            lowptjets = lowptjets[(lowptjets["pt_nomuon"] > 15)
                                  & (lowptjets["emef"] < 0.9)]
            # lowptjets lose their type here. Probably a bug, workaround
            lowptjets = ak.with_name(lowptjets, "PtEtaPhiMLorentzVector")
            if smear_met:
                lowptfac = lowptjets["juncfac"] * lowptjets["jerfac"] - 1
            else:
                lowptfac = lowptjets["juncfac"] - 1
            metx = metx - (
                ak.sum(jets.x * jets["factor"], axis=1)
                + ak.sum(lowptjets.x * lowptfac, axis=1))
            mety = mety - (
                ak.sum(jets.y * jets["factor"], axis=1)
                + ak.sum(lowptjets.y * lowptfac, axis=1))
        met["pt"] = np.hypot(metx, mety)
        met["phi"] = np.arctan2(mety, metx)
        return met

    def in_hem1516(self, phi, eta):
        """Return mask to select objects in faulty sector of hadronic
           calorimeter encap (HEM 15/16 issue in 2018 data)."""
        return ((-3.0 < eta) & (eta < -1.3) & (-1.57 < phi) & (phi < -0.87))

    def hem_cut(self, data):
        """Keep objects without the HEM 15/16 issue in 2018 data."""
        cut_ele = self.config["hem_cut_if_ele"]
        cut_muon = self.config["hem_cut_if_muon"]
        cut_jet = self.config["hem_cut_if_jet"]

        keep = np.full(len(data), True)
        if cut_ele:
            ele = data["Electron"]
            keep = keep & (~ak.any(self.in_hem1516(ele.phi, ele.eta)))
        if cut_muon:
            muon = data["Muon"]
            keep = keep & (~ak.any(self.in_hem1516(muon.phi, muon.eta)))
        if cut_jet:
            jet = data["Jet"]
            keep = keep & (~ak.any(self.in_hem1516(jet.phi, jet.eta)))
        return keep

    def lep_pt_requirement(self, data):
        """Require leptons with minimum pT threshold."""
        n = np.zeros(len(data))
        # This assumes leptons are ordered by pt highest first
        for i, pt_min in enumerate(self.config["lep_pt_min"]):
            # flatten to workaround awkward bug
            # https://github.com/scikit-hep/awkward-1.0/issues/1305
            mask = ak.flatten(ak.num(data["Lepton"]) > i, axis=None)
            n[mask] += np.asarray(
                pt_min < data["Lepton"].pt[mask][:, i]).astype(int)
        return n >= self.config["lep_pt_num_satisfied"]

    def good_mass_lepton_pair(self, data):
        """Which events have lepton pair mass required by the configuration"""
        return data["mll"] > self.config["mll_min"]

    def no_additional_leptons(self, is_mc, data):
        """Veto events with >= 3 leptons."""
        add_ele = self.electron_cuts(data["Electron"], good_lep=False)
        add_muon = self.muon_cuts(data["Muon"], good_lep=False)
        return ak.sum(add_ele, axis=1) + ak.sum(add_muon, axis=1) <= 2

    def has_jets(self, data):
        """Require events with minimum number of jets."""
        return self.config["num_jets_atleast"] <= ak.num(data["Jet"])

    def has_fatjets(self, data):
        """Require events with minimum number of fat jets."""
        return self.config["num_fatjets_atleast"] <= ak.num(data["FatJet"])

    def compute_puid_sys(self, central, weighter, wp, eta, pt,
                         pass_puid, has_gen_jet, sf_type):
        """Jet pileup ID systematic weights"""
        up = weighter(wp, eta, pt, pass_puid, has_gen_jet, sf_type, "up")
        down = weighter(wp, eta, pt, pass_puid, has_gen_jet, sf_type, "down")
        return (up / central, down / central)

    def jet_puid_sfs(self, data):
        """Compute event weights and systematics, if requested, for the jet
        pileup ID"""
        # Only apply SFs to jets for which the PU ID cut is applied,
        # i.e. pT < 50 GeV
        jets = data["Jet"][data["Jet"].pt < 50]
        wp = self.config["good_jet_puId"].split(":", 1)[1]
        pass_puid = jets["pass_pu_id"]
        has_gen_jet = jets["has_gen_jet"]
        eta = jets.eta
        pt = jets.pt
        weighter = self.config["jet_puid_sf"]
        systematics = {}
        weight = weighter(wp, eta, pt, pass_puid, has_gen_jet, "eff")
        if self.config["compute_systematics"]:
            systematics["jet_puid_eff"] = self.compute_puid_sys(
                weight, weighter, wp, eta, pt, pass_puid, has_gen_jet, "eff")
        if weighter.has_mis_prob:
            central = weighter(wp, eta, pt, pass_puid, has_gen_jet, "mis")
            weight = weight * central
            if self.config["compute_systematics"]:
                systematics["jet_puid_mis"] = self.compute_puid_sys(
                    weight, weighter, wp, eta, pt, pass_puid,
                    has_gen_jet, "mis")
        return weight, systematics

    def jet_pt_requirement(self, data):
        """Require jets with minimum pT threshold."""
        n = np.zeros(len(data))
        # This assumes jets are ordered by pt highest first
        for i, pt_min in enumerate(self.config["jet_pt_min"]):
            mask = ak.num(data["Jet"]) > i
            n[mask] += np.asarray(pt_min < data["Jet"].pt[mask][:, i]).astype(int)
        return n >= self.config["jet_pt_num_satisfied"]

    def fatjet_pt_requirement(self, data):
        """Require fat jets with minimum pT threshold."""
        n = np.zeros(len(data))
        # This assumes jets are ordered by pt highest first
        for i, pt_min in enumerate(self.config["fatjet_pt_min"]):
            mask = ak.num(data["FatJet"]) > i
            n[mask] += np.asarray(pt_min < data["FatJet"].pt[mask][:, i]).astype(int)
        return n >= self.config["fatjet_pt_num_satisfied"]

    def compute_btag_sys(self, central, up_name, down_name, weighter, wp, flav,
                         eta, pt, discr, efficiency):
        """b tagging systematic weights"""
        up = weighter(wp, flav, eta, pt, discr, up_name, efficiency)
        down = weighter(wp, flav, eta, pt, discr, down_name, efficiency)
        return (up / central, down / central)

    def compute_weight_btag(self, data, efficiency="central", never_sys=False):
        """Compute event weights and systematics, if requested, for the b
        tagging"""
        jets = data["Jet"]
        wp = self.config["btag"].split(":", 1)[1]
        flav = jets["hadronFlavour"]
        eta = jets.eta
        pt = jets.pt
        discr = jets["btag"]
        weight = np.ones(len(data))
        systematics = {}
        for i, weighter in enumerate(self.config["btag_sf"]):
            central = weighter(wp, flav, eta, pt, discr, "central", efficiency)
            if not never_sys and self.config["compute_systematics"]:
                if "btag_splitting_scheme" in self.config:
                    scheme = self.config["btag_splitting_scheme"].lower()
                elif ("split_btag_year_corr" in self.config and
                        self.config["split_btag_year_corr"]):
                    scheme = "years"
                else:
                    scheme = None
                if scheme is None:
                    light_unc_splits = heavy_unc_splits = {"": ""}
                elif scheme == "years":
                    light_unc_splits = heavy_unc_splits = \
                        {"corr": "_correlated", "uncorr": "_uncorrelated"}
                elif scheme == "sources":
                    heavy_unc_splits = {name: f"_{name}"
                                        for name in weighter.sources}
                    light_unc_splits = {"corr": "_correlated",
                                        "uncorr": "_uncorrelated"}
                else:
                    raise ValueError(
                        f"Invalid btag uncertainty scheme {scheme}")

                for name, split in heavy_unc_splits.items():
                    systematics[f"btagsf{i}" + name] = self.compute_btag_sys(
                        central, "heavy up" + split, "heavy down" + split,
                        weighter, wp, flav, eta, pt, discr, efficiency)
                for name, split in light_unc_splits.items():
                    systematics[f"btagsf{i}light" + name] = \
                        self.compute_btag_sys(
                            central, "light up" + split, "light down" + split,
                            weighter, wp, flav, eta, pt, discr, efficiency)
            weight = weight * central
        if never_sys:
            return weight
        else:
            return weight, systematics

    def scale_systematics_for_btag(self, selector, variation, dsname):
        """Modifies factors in the systematic table to account for differences
        in b-tag efficiencies. This is only done for variations the efficiency
        ROOT file contains a histogram with the name of the variation."""
        if "btag_sf" not in self.config or len(self.config["btag_sf"]) == 0:
            return
        available = set.intersection(
            *(w.available_efficiencies for w in self.config["btag_sf"]))
        data = selector.data
        systematics = selector.systematics
        central = self.compute_weight_btag(data, never_sys=True)
        if (variation == self.get_jetmet_nominal_arg()
                and dsname not in self.config["dataset_for_systematics"]):
            for name in ak.fields(systematics):
                if name == "weight":
                    continue
                if name not in available:
                    continue
                sys = systematics[name]
                varied_sf = self.compute_weight_btag(data, name, True)
                selector.set_systematic(name, sys / central * varied_sf)
        elif dsname in self.config["dataset_for_systematics"]:
            name = self.config["dataset_for_systematics"][dsname][1]
            if name in available:
                sys = systematics[name]
                varied_sf = self.compute_weight_btag(data, name, True)
                selector.set_systematic("weight", sys / central * varied_sf)
        elif variation.name in available:
            sys = systematics[name]
            varied_sf = self.compute_weight_btag(data, variation.name, True)
            selector.set_systematic("weight", sys / central * varied_sf)

    def btag_cut(self, is_mc, data):
        """Select events with minimum number of b-tagegd jets."""
        if self.config["btag"] == "skip":
            if self.config["num_atleast_btagged"] == 0:
                return np.full(len(data), True)
            else:
                raise RuntimeError(
                    "Cannot apply btag cut if btag is set to 'skip'. Please remove this "
                    "cut or choose a b-tagging algorithm.")
        num_btagged = ak.sum(data["Jet"]["btagged"], axis=1)
        accept = np.asarray(num_btagged >= self.config["num_atleast_btagged"])
        if is_mc and (
                "btag_sf" in self.config and len(self.config["btag_sf"]) != 0):
            weight, systematics = self.compute_weight_btag(data[accept])
            accept = accept.astype(float)
            accept[accept.astype(bool)] *= np.asarray(weight)
            return accept, systematics
        else:
            return accept

    def pick_lepton_pair(self, data):
        """Get one pair of leptons of opposite charge per event. The negative
        lepton comes first"""
        return data["Lepton"][
            ak.argsort(data["Lepton"]["pdgId"], ascending=False)]

    def pick_bs_from_lepton_pair(self, data):
        """Pick a bottom quark and a bottom antiquark that fit best to a pair
        of leptons assuming they come from a top pair decay. This is using
        the mlb histogram method"""
        recolepton = data["recolepton"]
        lep = recolepton[:, 0]
        antilep = recolepton[:, 1]
        # Build a reduced jet collection to avoid loading all branches and
        # make make this function faster overall
        columns = ["pt", "eta", "phi", "mass", "btagged"]
        jets = ak.with_name(data["Jet"][columns], "PtEtaPhiMLorentzVector")
        btags = jets[jets.btagged]
        jetsnob = jets[~jets.btagged]
        num_btags = ak.num(btags)
        b0, b1 = ak.unzip(ak.where(
            num_btags > 1, ak.combinations(btags, 2),
            ak.where(
                num_btags == 1, ak.cartesian([btags, jetsnob]),
                ak.combinations(jetsnob, 2))))
        bs = ak.concatenate([b0, b1], axis=1)
        bs_rev = ak.concatenate([b1, b0], axis=1)
        mass_alb = reduce(
            lambda a, b: a + b, ak.unzip(ak.cartesian([bs, antilep]))).mass
        mass_lb = reduce(
            lambda a, b: a + b, ak.unzip(ak.cartesian([bs_rev, lep]))).mass
        with uproot.open(self.config["reco_info_file"]) as f:
            mlb_prob = pepper.scale_factors.ScaleFactors.from_hist(f["mlb"])
        p_m_alb = mlb_prob(mlb=mass_alb)
        p_m_lb = mlb_prob(mlb=mass_lb)
        bestbpair_mlb = ak.unflatten(
            ak.argmax(p_m_alb * p_m_lb, axis=1), np.full(len(bs), 1))
        return ak.concatenate([bs[bestbpair_mlb], bs_rev[bestbpair_mlb]],
                              axis=1)

    def ttbar_system(self, reco_alg, rng, data):
        """Do ttbar reconstruction, obtaining four vectors for top pairs
        from their decay products"""
        lep = data["recolepton"][:, 0]
        antilep = data["recolepton"][:, 1]
        b = data["recob"][:, 0]
        antib = data["recob"][:, 1]
        met = data["MET"]
        if reco_alg == "sonnenschein":
            if self.config["reco_num_smear"] is None:
                energyfl = energyfj = 1
                alphal = alphaj = 0
                num_smear = 1
                mlb = None
            else:
                with uproot.open(self.config["reco_info_file"]) as f:
                    energyfl = f["energyfl"]
                    energyfj = f["energyfj"]
                    alphal = f["alphal"]
                    alphaj = f["alphaj"]
                    mlb = f["mlb"]
                num_smear = self.config["reco_num_smear"]
            if isinstance(self.config["reco_w_mass"], (int, float)):
                mw = self.config["reco_w_mass"]
            else:
                with uproot.open(self.config["reco_info_file"]) as f:
                    mw = f[self.config["reco_w_mass"]]
            if isinstance(self.config["reco_t_mass"], (int, float)):
                mt = self.config["reco_t_mass"]
            else:
                with uproot.open(self.config["reco_info_file"]) as f:
                    mt = f[self.config["reco_t_mass"]]
            top, antitop = sonnenschein(
                lep, antilep, b, antib, met, mwp=mw, mwm=mw, mt=mt, mat=mt,
                energyfl=energyfl, energyfj=energyfj, alphal=alphal,
                alphaj=alphaj, hist_mlb=mlb, num_smear=num_smear, rng=rng)
            return ak.concatenate([top, antitop], axis=1)
        else:
            raise ValueError(f"Invalid value for reco algorithm: {reco_alg}")

    def has_ttbar_system(self, data):
        """Whether the recontruction of the ttbar system was successful"""
        return ak.num(data["recot"]) > 0

    def build_nu_column_ttbar_system(self, data):
        """Get four momenta for the neutrinos coming from top pair decay"""
        lep = data["recolepton"][:, 0:1]
        antilep = data["recolepton"][:, 1:2]
        b = data["recob"][:, 0:1]
        antib = data["recob"][:, 1:2]
        top = data["recot"][:, 0:1]
        antitop = data["recot"][:, 1:2]
        nu = top - b - antilep
        antinu = antitop - antib - lep
        return ak.concatenate([nu, antinu], axis=1)
