import numpy as np
import awkward as ak


def veto_EEleak_flag(self, egammas: ak.Array) -> ak.Array:
    """
    Add branch to veto electrons/photons in the EE+ leak region for 2022. Ref to:
    https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis#Notes_on_addressing_EE_issue_in
    """
    if hasattr(egammas, "isScEtaEE"):
        # photons
        out_EEleak_region = (egammas.isScEtaEB) | (
            (egammas.isScEtaEE)
            & (egammas.seediEtaOriX < 45)
            & (egammas.seediPhiOriY > 72)
        )
    else:
        # electrons
        electron_ScEta = egammas.eta + egammas.deltaEtaSC
        out_EEleak_region = (np.abs(electron_ScEta) < self.gap_barrel_eta) | (
            (np.abs(electron_ScEta) > self.gap_endcap_eta)
            & (np.abs(electron_ScEta) < self.max_sc_eta)
            & (egammas.seediEtaOriX < 45)
            & (egammas.seediPhiOriY > 72)
        )

    egammas["vetoEELeak"] = out_EEleak_region

    return egammas
