import pepper
import awkward as ak
import numpy as np
from functools import partial
from coffea.nanoevents import PFNanoAODSchema


class Processor(pepper.ProcessorBasicPhysics):
    config_class = pepper.ConfigTTbarLL
    schema_class = PFNanoAODSchema

    def __init__(self, config, eventdir):
        # Need to call parent init to make histograms and such ready
        super().__init__(config, eventdir)

    def process_selection(self, selector, dsname, is_mc, filler):
        selector.set_multiple_columns(self.find_progenitor_quarks)
        self.unload_column("GenPart")
        selector.add_cut("Require_genparts", self.require_HP_genparts)
        selector.set_column("GenJet", self.calculate_genjet_pull)
        selector.set_multiple_columns(self.find_genjets_matching_HP)
        self.unload_column("GenJet")
        selector.add_cut("Require_matched_genjets", self.require_distinct_genjets)
        selector.set_column("Jet", self.calculate_jet_pull)
        selector.set_multiple_columns(self.find_jets_matching_HP)
        self.unload_column("Jet")
        selector.add_cut("Require_matched_jets", self.require_distinct_jets)
        selector.set_multiple_columns(self.set_pull_angles)

    def find_progenitor_quarks(self, data):
        new_cols = {}
        gen_parts = data["GenPart"]
        new_cols["gen_HP_b"] = ak.drop_none(
            gen_parts[(gen_parts.pdgId==5) & (gen_parts.parent.pdgId==6)]) # b quark from top quark
        new_cols["gen_HP_bbar"] = ak.drop_none(
            gen_parts[(gen_parts.pdgId==-5) & (gen_parts.parent.pdgId==-6)]) # Anti-b from antitop
        new_cols["gen_HP_qfromWplus"] = ak.drop_none(
            gen_parts[(abs(gen_parts.pdgId)<6) & (gen_parts.parent.pdgId==24) & (gen_parts.parent.distinctParent.pdgId==6)]) # q and qbar from W+ (originally from top)
        new_cols["gen_HP_qfromWminus"] = ak.drop_none(
            gen_parts[(abs(gen_parts.pdgId)<6) & (gen_parts.parent.pdgId==-24) & (gen_parts.parent.distinctParent.pdgId==-6)]) # q and qbar from W- (originally from antitop)
        return new_cols

    def require_HP_genparts(self, data):
        return ((ak.num(data["gen_HP_b"]) == 1) &
                (ak.num(data["gen_HP_bbar"]) == 1) &
                (ak.num(data["gen_HP_qfromWplus"]) == 2) &
                (ak.num(data["gen_HP_qfromWminus"]) == 2))

    def calculate_genjet_pull(self, data):
        def weighted_sum(deltas, pt_cands, pt_jet, mag):
            return ak.sum(deltas * pt_cands * mag / pt_jet, axis=2)

        jets = data["GenJet"]
        pairs = ak.cartesian({"jet": jets, "gp": data["GenCands"]}, nested=True)
        dr = pairs.jet.delta_r(pairs.gp)
        mask = dr < 0.4
        cands = pairs.gp[mask]
        mag = cands.deltaRapidityPhi(jets)
        dphis = cands.deltaphi(jets)
        dy = cands.rapidity - jets.rapidity
        jets["pull_phi"] = weighted_sum(dphis, cands.pt, jets.pt, mag)
        jets["pull_rapidity"] = weighted_sum(dy, cands.pt, jets.pt, mag)
        # Some jets don't have matched Gen candidates, like due to the choices made when making this nano - drop these for now
        jets = jets[ak.num(cands, axis=2) > 0]
        return jets

    def calculate_jet_pull(self, data):
        def weighted_sum(deltas, pt_cands, pt_jet, mag):
            return ak.sum(deltas * pt_cands * mag / pt_jet, axis=2)

        jets = data["Jet"]
        #cands = jets.candidates
        cands = jets.constituents.pf
        mag = cands.deltaRapidityPhi(jets)
        dphis = cands.deltaphi(jets)
        dy = cands.rapidity - jets.rapidity
        jets["pull_phi"] = weighted_sum(dphis, cands.pt, jets.pt, mag)
        jets["pull_rapidity"] = weighted_sum(dy, cands.pt, jets.pt, mag)
        jets = jets[ak.num(cands, axis=2) > 0]
        return jets

    def find_genjets_matching_HP(self, data):
        new_cols = {}
        gen_jets = data["GenJet"]
        gen_jets = ak.with_field(gen_jets, ak.local_index(gen_jets, axis=1), "idx")
        new_cols["genjet_from_HP_b"] = ak.drop_none(data["gen_HP_b"].nearest(gen_jets, threshold=0.3))
        new_cols["genjet_from_HP_bbar"] = ak.drop_none(data["gen_HP_bbar"].nearest(gen_jets, threshold=0.3))
        new_cols["genjet_from_HP_qfromWplus"] = ak.drop_none(data["gen_HP_qfromWplus"].nearest(gen_jets, threshold=0.3))
        new_cols["genjet_from_HP_qfromWminus"] = ak.drop_none(data["gen_HP_qfromWminus"].nearest(gen_jets, threshold=0.3))
        return new_cols

    def find_jets_matching_HP(self, data):
        new_cols = {}
        jets = data["Jet"]
        jets = ak.with_field(jets, ak.local_index(jets, axis=1), "idx")
        new_cols["jet_from_HP_b"] = ak.drop_none(data["gen_HP_b"].nearest(jets, threshold=0.3))
        new_cols["jet_from_HP_bbar"] = ak.drop_none(data["gen_HP_bbar"].nearest(jets, threshold=0.3))
        new_cols["jet_from_HP_qfromWplus"] = ak.drop_none(data["gen_HP_qfromWplus"].nearest(jets, threshold=0.3))
        new_cols["jet_from_HP_qfromWminus"] = ak.drop_none(data["gen_HP_qfromWminus"].nearest(jets, threshold=0.3))
        return new_cols

    def require_distinct_genjets(self, data):
        idx_all = ak.concatenate(
            [j.idx for j in [data["genjet_from_HP_b"], data["genjet_from_HP_bbar"], data["genjet_from_HP_qfromWplus"], data["genjet_from_HP_qfromWminus"]]], axis=1)
        idx_sorted = ak.sort(idx_all, axis=1)
        all_distinct = ak.all(idx_sorted[:, 1:] != idx_sorted[:, :-1], axis=1)
        return ((ak.num(data["genjet_from_HP_b"]) == 1) &
                (ak.num(data["genjet_from_HP_bbar"]) == 1) &
                (ak.num(data["genjet_from_HP_qfromWplus"]) == 2) &
                (ak.num(data["genjet_from_HP_qfromWminus"]) == 2) &
                all_distinct)

    def require_distinct_jets(self, data):
        idx_all = ak.concatenate(
            [j.idx for j in [data["jet_from_HP_b"], data["jet_from_HP_bbar"], data["jet_from_HP_qfromWplus"], data["jet_from_HP_qfromWminus"]]], axis=1)
        idx_sorted = ak.sort(idx_all, axis=1)
        all_distinct = ak.all(idx_sorted[:, 1:] != idx_sorted[:, :-1], axis=1)
        return ((ak.num(data["jet_from_HP_b"]) == 1) &
                (ak.num(data["jet_from_HP_bbar"]) == 1) &
                (ak.num(data["jet_from_HP_qfromWplus"]) == 2) &
                (ak.num(data["jet_from_HP_qfromWminus"]) == 2) &
                all_distinct)

    @staticmethod
    def rectify_angle(phi):
        return (phi + np.pi) % (2 * np.pi) - np.pi

    def calculate_angle(self, phi1, y1, phi2, y2):
        return self.rectify_angle(np.arctan2(phi1*y2 - phi2*y1, phi1*phi2 + y1*y2))

    def calculate_pull_angle(self, jets):
        jcv_y = jets[:, 1].rapidity -  jets[:, 0].rapidity # JCV = "jet connection vector"
        jcv_phi = jets[:, 1].delta_phi(jets[:, 0])
        forward_pull_angle = self.calculate_angle(jcv_phi, jcv_y, jets["pull_phi"][:, 0], jets["pull_rapidity"][:, 0])
        backward_pull_angle = self.calculate_angle(-jcv_phi, -jcv_y, jets["pull_phi"][:, 1], jets["pull_rapidity"][:, 1])
        return forward_pull_angle, backward_pull_angle

    def set_pull_angles(self, data):
        new_cols = {}
        new_cols["forward_pull_angle_Wplus"], new_cols["backward_pull_angle_Wplus"] = \
            self.calculate_pull_angle(data["genjet_from_HP_qfromWplus"])
        new_cols["forward_pull_angle_Wminus"], new_cols["backward_pull_angle_Wminus"] = \
            self.calculate_pull_angle(data["genjet_from_HP_qfromWminus"])
        new_cols["forward_pull_angle_bs"], new_cols["backward_pull_angle_bs"] = \
            self.calculate_pull_angle(ak.concatenate([data["genjet_from_HP_b"], data["genjet_from_HP_bbar"]], axis=1))
        new_cols["reco_forward_pull_angle_Wplus"], new_cols["reco_backward_pull_angle_Wplus"] = \
            self.calculate_pull_angle(data["jet_from_HP_qfromWplus"])
        new_cols["reco_forward_pull_angle_Wminus"], new_cols["reco_backward_pull_angle_Wminus"] = \
            self.calculate_pull_angle(data["jet_from_HP_qfromWminus"])
        new_cols["reco_forward_pull_angle_bs"], new_cols["reco_backward_pull_angle_bs"] = \
            self.calculate_pull_angle(ak.concatenate([data["jet_from_HP_b"], data["jet_from_HP_bbar"]], axis=1))
        return new_cols


