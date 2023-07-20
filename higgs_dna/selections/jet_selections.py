from higgs_dna.selections.object_selections import delta_r_mask
import awkward


def select_jets(
    self,
    jets: awkward.highlevel.Array,
    diphotons: awkward.highlevel.Array,
    muons: awkward.highlevel.Array,
    electrons: awkward.highlevel.Array,
) -> awkward.highlevel.Array:
    pt_cut = jets.pt > self.jet_pt_threshold
    eta_cut = abs(jets.eta) < self.jet_max_eta
    dr_dipho_cut = awkward.ones_like(pt_cut) > 0
    if self.clean_jet_dipho & (awkward.count(diphotons) > 0):
        dr_dipho_cut = delta_r_mask(jets, diphotons, self.jet_dipho_min_dr)

    if (self.clean_jet_pho) & (awkward.count(diphotons) > 0):
        lead = awkward.zip(
            {
                "pt": diphotons.pho_lead.pt,
                "eta": diphotons.pho_lead.eta,
                "phi": diphotons.pho_lead.phi,
                "mass": diphotons.pho_lead.mass,
                "charge": diphotons.pho_lead.charge,
            }
        )
        lead = awkward.with_name(lead, "PtEtaPhiMCandidate")
        sublead = awkward.zip(
            {
                "pt": diphotons.pho_sublead.pt,
                "eta": diphotons.pho_sublead.eta,
                "phi": diphotons.pho_sublead.phi,
                "mass": diphotons.pho_sublead.mass,
                "charge": diphotons.pho_sublead.charge,
            }
        )
        sublead = awkward.with_name(sublead, "PtEtaPhiMCandidate")
        dr_pho_lead_cut = delta_r_mask(jets, lead, self.jet_pho_min_dr)
        dr_pho_sublead_cut = delta_r_mask(jets, sublead, self.jet_pho_min_dr)
    else:
        dr_pho_lead_cut = jets.pt > -1
        dr_pho_sublead_cut = jets.pt > -1

    if (self.clean_jet_ele) & (awkward.count(electrons) > 0):
        dr_electrons_cut = delta_r_mask(jets, electrons, self.jet_ele_min_dr)
    else:
        dr_electrons_cut = jets.pt > -1

    if (self.clean_jet_muo) & (awkward.count(muons) > 0):
        dr_muons_cut = delta_r_mask(jets, muons, self.jet_muo_min_dr)
    else:
        dr_muons_cut = jets.pt > -1

    return (
        (pt_cut)
        & (eta_cut)
        & (dr_dipho_cut)
        & (dr_pho_lead_cut)
        & (dr_pho_sublead_cut)
        & (dr_electrons_cut)
        & (dr_muons_cut)
    )
