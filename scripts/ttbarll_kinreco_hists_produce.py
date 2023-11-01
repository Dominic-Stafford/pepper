#!/usr/bin/env python3

import os
import numpy as np
import awkward as ak
import uproot

import pepper


class Processor(pepper.ProcessorTTbarLL):
    def __init__(self, config, destdir):
        if "reco_algorithm" in config:
            del config["reco_algorithm"]
        if "reco_info_file" in config:
            del config["reco_info_file"]
        if "blinding_denom" in config:
            del config["blinding_denom"]
        config["compute_systematics"] = False
        config["histogram_format"] = "hist"
        config["hists"] = {
            "mlb": self._make_hist(
                "mlb", r"$m_{\mathrm{lb}}$ (GeV)", 200, 0, 200, "mlb"),
            "mw": self._make_hist(
                "mw", r"$m_{\mathrm{W}}$ (GeV)", 160, 40, 120, "mw"),
            "mt": self._make_hist(
                "mt", r"$m_{\mathrm{t}}$ (GeV)", 60, 160, 190, "mt"),
            "alphal": self._make_hist(
                "alpha", r"$\alpha$ (rad)", 20, 0, 0.02, "alphal"),
            "alphaj": self._make_hist(
                "alpha", r"$\alpha$ (rad)", 100, 0, 0.2, "alphaj"),
            "energyfl": self._make_hist(
                "energyf", r"$E_{\mathrm{gen}} / E_{\mathrm{reco}}$", 200, 0.5,
                1.5, "energyfl"),
            "energyfj": self._make_hist(
                "energyf", r"$E_{\mathrm{gen}} / E_{\mathrm{reco}}$", 200, 0,
                3, "energyfj"),
        }

        super().__init__(config, None)

    @staticmethod
    def _make_hist(name, label, n_or_arr, lo, hi, fill):
        return pepper.HistDefinition({
            "bins": [
                {
                    "name": name,
                    "label": label,
                    "n_or_arr": n_or_arr,
                    "lo": lo,
                    "hi": hi
                }
            ],
            "fill": {
                name: [
                    fill
                ]
            }
        })

    def preprocess(self, datasets):
        return {"TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8":
                datasets["TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8"]}

    def process_selection(self, selector, dsname, is_mc, filler):
        selector.set_multiple_columns(self.build_gen_columns)
        selector.add_cut("HasGenParticles", self.has_gen_particles)
        selector.set_column("mlb", self.mlb)
        selector.set_column("mw", self.mw)
        selector.set_column("mt", self.mt)

        super().process_selection(selector, dsname, is_mc, filler)

        selector.set_multiple_columns(self.lepton_recogen)
        selector.set_multiple_columns(self.jet_recogen)

    @staticmethod
    def sortby(data, field):
        sorting = ak.argsort(data[field], ascending=False)
        return data[sorting]

    def build_gen_columns(self, data):
        part = data["GenPart"]
        mass = np.asarray(ak.flatten(part.mass))
        # mass is only stored for masses greater than 10 GeV
        # this is a fudge to set the b mass to the right
        # value if not stored
        mass[(mass == 0) & (abs(ak.flatten(part["pdgId"])) == 5)] = 4.18
        mass[(mass == 0) & (abs(ak.flatten(part["pdgId"])) == 13)] = 0.102
        part["mass"] = ak.unflatten(mass, ak.num(part))
        # Fill none in order to not enter masked arrays regime
        part["motherid"] = ak.fill_none(part.parent.pdgId, 0)
        part = part[part["motherid"] != 0]
        part = part[part.hasFlags("isFirstCopy", "isHardProcess")]
        abspdg = abs(part.pdgId)
        sgn = np.sign(part.pdgId)

        cols = {}
        cols["genlepton"] = part[((abspdg == 11) | (abspdg == 13))
                                 & (part.motherid == sgn * -24)]
        cols["genlepton"] = self.sortby(cols["genlepton"], "pdgId")

        cols["genb"] = part[(abspdg == 5) & (part.motherid == sgn * 6)]
        cols["genb"] = self.sortby(cols["genb"], "pdgId")

        cols["genw"] = part[(abspdg == 24) & (part.motherid == sgn * 6)]
        cols["genw"] = self.sortby(cols["genw"], "pdgId")

        cols["gent"] = part[(abspdg == 6)]
        cols["gent"] = self.sortby(cols["gent"], "pdgId")

        return cols

    def has_gen_particles(self, data):
        return ((ak.num(data["genlepton"]) == 2)
                & (ak.num(data["genb"]) == 2)
                & (ak.num(data["genw"]) == 2)
                & (ak.num(data["gent"]) == 2))

    def mlb(self, data):
        return (data["genlepton"][:, ::-1] + data["genb"]).mass

    def mw(self, data):
        return data["genw"].mass

    def mt(self, data):
        return data["gent"].mass

    def match_leptons(self, data):
        recolep = self.sortby(data["Lepton"][:, :2], "pdgId")
        genlep = data["genlepton"]
        is_same_flavor = recolep.pdgId == genlep.pdgId
        is_close = recolep.delta_r(genlep) < 0.3
        is_matched = is_same_flavor & is_close

        return genlep[is_matched], recolep[is_matched]

    def match_jets(self, data):
        recojet = ak.with_name(data["Jet"][["pt", "eta", "phi", "mass"]],
                               "PtEtaPhiMLorentzVector")
        genb = ak.with_name(data["genb"][["pt", "eta", "phi", "mass"]],
                            "PtEtaPhiMLorentzVector")
        genbc, recojetc = ak.unzip(ak.cartesian([genb, recojet], nested=True))
        is_close = genbc.delta_r(recojetc) < 0.3
        is_matched = is_close & (ak.sum(is_close, axis=2) == 1)

        mrecojet = [recojet[is_matched[:, i]] for i in range(2)]
        mrecojet = ak.concatenate(mrecojet, axis=1)
        return genb[ak.any(is_matched, axis=2)], mrecojet

    def lepton_recogen(self, data):
        gen, reco = self.match_leptons(data)
        energyf = gen.energy / reco.energy
        deltaphi = gen.delta_phi(reco)
        return {"energyfl": energyf, "alphal": deltaphi}

    def jet_recogen(self, data):
        gen, reco = self.match_jets(data)
        energyf = gen.energy / reco.energy
        deltaphi = gen.delta_phi(reco)
        return {"energyfj": energyf, "alphaj": deltaphi}

    def save_output(self, output, dest):
        output = output["hists"]["TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8"]
        items = ("mlb", "mw", "mt", "alphal", "energyfl", "alphaj",
                 "energyfj")
        cuts = ("HasGenParticles", "ReqMET")
        with uproot.recreate(os.path.join(dest, "kinreco.root")) as f:
            for item in items:
                for cut in cuts:
                    if (cut, item) in output:
                        axname = item[:-1] if item[-1] in "jl" else item
                        f[item] = output[(cut, item)].project(axname)
                        break
                else:
                    raise RuntimeError(f"No histogram for {item}")


if __name__ == "__main__":
    from pepper import runproc
    runproc.run_processor(Processor, mconly=True)
