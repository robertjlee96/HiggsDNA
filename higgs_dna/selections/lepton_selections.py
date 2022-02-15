from higgs_dna.selections.object_selections import delta_r_mask
import awkward  # type: ignore


def select_electrons(
    electrons: awkward.highlevel.Array,
    diphotons: awkward.highlevel.Array,
    electron_pt_threshold: float,
) -> awkward.highlevel.Array:
    pt_cut = electrons.pt > electron_pt_threshold
    eta_cut = abs(electrons.eta) < 2.4
    id_cut = electrons.mvaFall17V2Iso_WP90
    dr_pho_cut = delta_r_mask(electrons, diphotons, 0.2)

    return pt_cut & eta_cut & id_cut & dr_pho_cut


def select_muons(
    muons: awkward.highlevel.Array, diphotons: awkward.highlevel.Array
) -> awkward.highlevel.Array:
    pt_cut = muons.pt > 25
    eta_cut = abs(muons.eta) < 2.4
    id_cut = muons.mediumId
    dr_pho_cut = delta_r_mask(muons, diphotons, 0.2)

    return pt_cut & eta_cut & id_cut & dr_pho_cut
