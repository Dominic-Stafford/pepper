import pepper
import awkward as ak
import numpy as np
import fastjet
import vector
from functools import partial
from coffea.nanoevents import PFNanoAODSchema

vector.register_awkward()


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
        selector.set_column("genjet_connected_idx", self.assign_connected_genjets)
        selector.set_multiple_columns(self.count_correct_connections)
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
        # Same pull again, but measured from the winner-take-all axis
        wta_y, wta_phi = self.calculate_wta_axis(cands)
        dphis_wta = self.rectify_angle(cands.phi - wta_phi)
        dy_wta = cands.rapidity - wta_y
        mag_wta = np.sqrt(dphis_wta**2 + dy_wta**2)
        jets["wta_rapidity"] = wta_y
        jets["wta_phi"] = wta_phi
        jets["pull_phi_wta"] = weighted_sum(dphis_wta, cands.pt, jets.pt, mag_wta)
        jets["pull_rapidity_wta"] = weighted_sum(dy_wta, cands.pt, jets.pt, mag_wta)
        # Some jets don't have matched Gen candidates, like due to the choices made when making this nano - drop these for now
        jets = jets[ak.num(cands, axis=2) > 0]
        return jets

    # Reclustering definition for the winner-take-all axis. R is larger than
    # the AK4 jets so that all constituents end up in a single jet.
    wta_jetdef = fastjet.JetDefinition(
        fastjet.antikt_algorithm, 0.8, fastjet.WTA_pt_scheme)

    def calculate_wta_axis(self, cands):
        """Recluster the constituents of each jet with the winner-take-all
        recombination scheme and return the (rapidity, phi) of the resulting axis"""
        njets = ak.num(cands, axis=1)
        flat = ak.flatten(cands, axis=1)
        inp = ak.zip({"pt": flat.pt, "eta": flat.eta, "phi": flat.phi,
                      "mass": flat.mass}, with_name="Momentum4D")
        subjets = fastjet.ClusterSequence(inp, self.wta_jetdef).inclusive_jets()
        subjets = ak.with_name(subjets, "Momentum4D")
        # Normally everything ends up in one subjet - take the hardest if not
        hardest = ak.firsts(
            subjets[ak.argsort(subjets.pt, axis=1, ascending=False)], axis=1)
        return (ak.unflatten(hardest.rapidity, njets, axis=0),
                ak.unflatten(hardest.phi, njets, axis=0))

    def calculate_jet_pull(self, data):
        def weighted_sum(deltas, pt_cands, pt_jet, mag):
            return ak.sum(deltas * pt_cands * mag / pt_jet, axis=2)

        jets = data["Jet"]
        cands = jets.constituents.pf
        mag = cands.deltaRapidityPhi(jets)
        dphis = cands.deltaphi(jets)
        dy = cands.rapidity - jets.rapidity
        jets["pull_phi"] = weighted_sum(dphis, cands.pt, jets.pt, mag)
        jets["pull_rapidity"] = weighted_sum(dy, cands.pt, jets.pt, mag)
        # Same pull, but each candidate scaled by its PUPPI weight.
        # The denominator is the PUPPI-weighted pt sum so the weights still add up to 1.
        pt_puppi = cands.pt * cands.puppiWeight
        pt_jet_puppi = ak.sum(pt_puppi, axis=2)
        jets["pull_phi_puppi"] = weighted_sum(dphis, pt_puppi, pt_jet_puppi, mag)
        jets["pull_rapidity_puppi"] = weighted_sum(dy, pt_puppi, pt_jet_puppi, mag)
        # Same pull again, but measured from the winner-take-all axis instead of the jet axis.
        wta_y, wta_phi = self.calculate_wta_axis(cands)
        dphis_wta = self.rectify_angle(cands.phi - wta_phi)
        dy_wta = cands.rapidity - wta_y
        mag_wta = np.sqrt(dphis_wta**2 + dy_wta**2)
        jets["wta_rapidity"] = wta_y
        jets["wta_phi"] = wta_phi
        jets["pull_phi_wta"] = weighted_sum(dphis_wta, cands.pt, jets.pt, mag_wta)
        jets["pull_rapidity_wta"] = weighted_sum(dy_wta, cands.pt, jets.pt, mag_wta)
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

    @staticmethod
    def get_W_genjets(data):
        """The four W decay jets, ordered [W+, W+, W-, W-], so that the true
        colour connected partner of a jet is the other one of its pair."""
        return ak.concatenate([data["genjet_from_HP_qfromWplus"],
                               data["genjet_from_HP_qfromWminus"]], axis=1)

    def assign_connected_genjets(self, data):
        """For each of the four W decay jets, find the one its pull vector
        points at most closely. Returns the position (0 to 3) of that jet."""
        jets = self.get_W_genjets(data)
        jets = ak.with_field(jets, ak.local_index(jets, axis=1), "pos")
        pairs = ak.cartesian({"i": jets, "j": jets}, nested=True)
        jcv_y = pairs.j.rapidity - pairs.i.rapidity
        jcv_phi = self.rectify_angle(pairs.j.phi - pairs.i.phi)
        angle = self.calculate_angle(jcv_phi, jcv_y, pairs.i["pull_phi"], pairs.i["pull_rapidity"])
        # Mask out pairing a jet with itself, then take the smallest angle
        angle = ak.mask(abs(angle), pairs.i.pos != pairs.j.pos)
        closest = ak.argmin(angle, axis=2, keepdims=True)
        return ak.firsts(pairs.j.pos[closest], axis=2)

    def count_correct_connections(self, data):
        """Count how many of the four W decay jets got assigned their true
        colour connected partner, i.e. the other jet from the same W."""
        connected = data["genjet_connected_idx"]
        pos = ak.local_index(connected, axis=1)
        true_partner = ak.where(pos % 2 == 0, pos + 1, pos - 1)
        return {"n_correct_connections": ak.sum(connected == true_partner, axis=1)}

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

    def calculate_pull_angle(self, jets, phi_field="pull_phi", y_field="pull_rapidity",
                             axis_y_field=None, axis_phi_field=None):
        if axis_y_field is None:
            jcv_y = jets[:, 1].rapidity -  jets[:, 0].rapidity # JCV = "jet connection vector"
            jcv_phi = jets[:, 1].delta_phi(jets[:, 0])
        else:
            # Build the jet connection vector from a different axis, e.g. WTA
            jcv_y = jets[axis_y_field][:, 1] - jets[axis_y_field][:, 0]
            jcv_phi = self.rectify_angle(
                jets[axis_phi_field][:, 1] - jets[axis_phi_field][:, 0])
        forward_pull_angle = self.calculate_angle(jcv_phi, jcv_y, jets[phi_field][:, 0], jets[y_field][:, 0])
        backward_pull_angle = self.calculate_angle(-jcv_phi, -jcv_y, jets[phi_field][:, 1], jets[y_field][:, 1])
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
        puppi = {"phi_field": "pull_phi_puppi", "y_field": "pull_rapidity_puppi"}
        new_cols["puppi_forward_pull_angle_Wplus"], new_cols["puppi_backward_pull_angle_Wplus"] = \
            self.calculate_pull_angle(data["jet_from_HP_qfromWplus"], **puppi)
        new_cols["puppi_forward_pull_angle_Wminus"], new_cols["puppi_backward_pull_angle_Wminus"] = \
            self.calculate_pull_angle(data["jet_from_HP_qfromWminus"], **puppi)
        new_cols["puppi_forward_pull_angle_bs"], new_cols["puppi_backward_pull_angle_bs"] = \
            self.calculate_pull_angle(ak.concatenate([data["jet_from_HP_b"], data["jet_from_HP_bbar"]], axis=1), **puppi)
        wta = {"phi_field": "pull_phi_wta", "y_field": "pull_rapidity_wta",
               "axis_y_field": "wta_rapidity", "axis_phi_field": "wta_phi"}
        new_cols["genwta_forward_pull_angle_Wplus"], new_cols["genwta_backward_pull_angle_Wplus"] = \
            self.calculate_pull_angle(data["genjet_from_HP_qfromWplus"], **wta)
        new_cols["genwta_forward_pull_angle_Wminus"], new_cols["genwta_backward_pull_angle_Wminus"] = \
            self.calculate_pull_angle(data["genjet_from_HP_qfromWminus"], **wta)
        new_cols["genwta_forward_pull_angle_bs"], new_cols["genwta_backward_pull_angle_bs"] = \
            self.calculate_pull_angle(ak.concatenate([data["genjet_from_HP_b"], data["genjet_from_HP_bbar"]], axis=1), **wta)
        new_cols["wta_forward_pull_angle_Wplus"], new_cols["wta_backward_pull_angle_Wplus"] = \
            self.calculate_pull_angle(data["jet_from_HP_qfromWplus"], **wta)
        new_cols["wta_forward_pull_angle_Wminus"], new_cols["wta_backward_pull_angle_Wminus"] = \
            self.calculate_pull_angle(data["jet_from_HP_qfromWminus"], **wta)
        new_cols["wta_forward_pull_angle_bs"], new_cols["wta_backward_pull_angle_bs"] = \
            self.calculate_pull_angle(ak.concatenate([data["jet_from_HP_b"], data["jet_from_HP_bbar"]], axis=1), **wta)
        return new_cols
