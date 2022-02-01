from higgs_dna.selections.object_selections import delta_r_mask
import awkward as ak


def select_jets(jets, diphotons, muons, electrons):
    pt_cut = jets.pt > 25
    eta_cut = abs(jets.eta) < 2.4
    dr_pho_cut = delta_r_mask(jets, diphotons, 0.4)
    dr_electrons_cut = delta_r_mask(jets, electrons, 0.4)
    dr_muons_cut = delta_r_mask(jets, muons, 0.4)

    return pt_cut & eta_cut & dr_pho_cut & dr_electrons_cut & dr_muons_cut
