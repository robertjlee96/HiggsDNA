from higgs_dna.selections.object_selections import delta_r_mask
import awkward


def select_electrons(
    self,
    electrons: awkward.highlevel.Array,
    diphotons: awkward.highlevel.Array,
) -> awkward.highlevel.Array:
    pt_cut = electrons.pt > self.electron_pt_threshold

    eta_cut = abs(electrons.eta) < self.electron_max_eta

    if self.el_iso_wp == "WP90":
        id_cut = electrons.mvaIso_WP90
    elif self.el_iso_wp == "WP80":
        id_cut = electrons.mvaIso_WP80
    # WPL is not supported anymore with the Run 3 electron ID, CMSSW 130X v12 nanoAODs only have WP80 and WP90 options
    # elif self.el_iso_wp == "WPL":
    #    id_cut = electrons.mvaIso_WPL
    else:
        id_cut = electrons.pt > 0.

    dr_pho_cut = delta_r_mask(electrons, diphotons, 0.2)

    return pt_cut & eta_cut & id_cut & dr_pho_cut


def select_muons(
    self,
    muons: awkward.highlevel.Array,
    diphotons: awkward.highlevel.Array
) -> awkward.highlevel.Array:
    pt_cut = muons.pt > self.muon_pt_threshold

    eta_cut = abs(muons.eta) < self.muon_max_eta

    if self.mu_iso_wp == "tight":
        id_cut = muons.tightId
    elif self.mu_iso_wp == "medium":
        id_cut = muons.mediumId
    elif self.mu_iso_wp == "loose":
        id_cut = muons.looseId
    else:
        id_cut = muons.pt > 0

    if self.global_muon:
        global_cut = muons.isGlobal
    else:
        global_cut = muons.pt > 0

    dr_pho_cut = delta_r_mask(muons, diphotons, 0.2)

    return pt_cut & eta_cut & id_cut & dr_pho_cut & global_cut
