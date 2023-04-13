import awkward
import numpy


# photon preselection -> take as input nAOD Photon collection and return the Photons that pass
# cuts (pt, eta, sieie, mvaID, iso... etc)
#
def photon_preselection(
    self, photons: awkward.Array, events: awkward.Array
) -> awkward.Array:
    """
    Apply preselection cuts to photons.
    Note that these selections are applied on each photon, it is not based on the diphoton pair.
    """
    # hlt-mimicking cuts
    rho = events.Rho.fixedGridRhoAll * awkward.ones_like(photons.pt)
    photon_abs_eta = numpy.abs(photons.eta)
    isEB_high_r9 = (photon_abs_eta < self.gap_barrel_eta) & (
        photons.r9 > self.min_full5x5_r9_EB_high_r9
    )
    isEE_high_r9 = (photon_abs_eta > self.gap_endcap_eta) & (
        photons.r9 > self.min_full5x5_r9_EE_high_r9
    )
    isEB_low_r9 = (
        (photon_abs_eta < self.gap_barrel_eta)
        & (photons.r9 > self.min_full5x5_r9_EB_low_r9)
        & (photons.r9 < self.min_full5x5_r9_EB_high_r9)
        & (
            photons.pfChargedIsoPFPV  # for v11
            < self.max_trkSumPtHollowConeDR03_EB_low_r9
        )
        & (photons.sieie < self.max_sieie_EB_low_r9)
        & (
            (
                (photon_abs_eta < self.eta_rho_corr)
                & (
                    photons.pfPhoIso03 - rho * self.low_eta_rho_corr
                    < self.max_pho_iso_EB_low_r9
                )
            )
            | (
                (photon_abs_eta > self.eta_rho_corr)
                & (
                    photons.pfPhoIso03 - rho * self.high_eta_rho_corr
                    < self.max_pho_iso_EB_low_r9
                )
            )
        )
    )
    isEE_low_r9 = (
        (photon_abs_eta < self.gap_barrel_eta)
        & (photons.r9 > self.min_full5x5_r9_EE_low_r9)
        & (photons.r9 < self.min_full5x5_r9_EE_high_r9)
        & (
            photons.pfChargedIsoPFPV  # for v11
            < self.max_trkSumPtHollowConeDR03_EE_low_r9
        )
        & (photons.sieie < self.max_sieie_EE_low_r9)
        & (
            (
                (photon_abs_eta < self.eta_rho_corr)
                & (
                    photons.pfPhoIso03 - rho * self.low_eta_rho_corr
                    < self.max_pho_iso_EE_low_r9
                )
            )
            | (
                (photon_abs_eta > self.eta_rho_corr)
                & (
                    photons.pfPhoIso03 - rho * self.high_eta_rho_corr
                    < self.max_pho_iso_EE_low_r9
                )
            )
        )
    )

    return photons[
        (photons.electronVeto > self.e_veto)
        & (photons.pt > self.min_pt_photon)
        & (photon_abs_eta < self.max_sc_eta)
        & (
            (photon_abs_eta < self.gap_barrel_eta)
            | (photon_abs_eta > self.gap_endcap_eta)
        )
        & (photons.mvaID > self.min_mvaid)
        & (photons.hoe < self.max_hovere)
        & (
            (photons.r9 > self.min_full5x5_r9)
            | (
                photons.pfRelIso03_chg_quadratic < self.max_chad_iso
            )  # changed from pfRelIso03_chg since this variable is not in v11 nanoAOD...?
            | (photons.pfRelIso03_chg_quadratic / photons.pt < self.max_chad_rel_iso)
        )
        & (isEB_high_r9 | isEB_low_r9 | isEE_high_r9 | isEE_low_r9)
    ]
