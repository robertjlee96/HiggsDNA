from higgs_dna.utils import awkward_utils
from higgs_dna.systematics.utils import systematic_from_bins, ic_systematic_from_bins
from higgs_dna.systematics.data.electron_veto_sf import (
    PHOTON_ELECTRON_VETO_SF_2016,
    PHOTON_ELECTRON_VETO_SF_2017,
    PHOTON_ELECTRON_VETO_SF_2018,
)
from higgs_dna.systematics.data.trigger_sf import (
    LEAD_TRIGGER_SF_2016,
    SUBLEAD_TRIGGER_SF_2016,
    LEAD_TRIGGER_SF_2017,
    SUBLEAD_TRIGGER_SF_2017,
    LEAD_TRIGGER_SF_2018,
    SUBLEAD_TRIGGER_SF_2018,
)
from higgs_dna.systematics.data.material import (
    MATERIAL_2016,
    MATERIAL_2017,
    MATERIAL_2018,
)

from higgs_dna.systematics.data.fnuf import FNUF_2016, FNUF_2017, FNUF_2018

import awkward
import logging

logger = logging.getLogger(__name__)

PHOTON_PRESELECTION_SFs = {
    "variables": ["photon_eta", "photon_r9"],
    "bins": [
        {
            "photon_eta": [0.0, 1.5],
            "photon_r9": [0.0, 0.85],
            "value": 1.0057,
            "uncertainty": 0.0010,
        },
        {
            "photon_eta": [0.0, 1.5],
            "photon_r9": [0.85, 999.0],
            "value": 0.9988,
            "uncertainty": 0.0009,
        },
        {
            "photon_eta": [1.5, 6.0],
            "photon_r9": [0.0, 0.9],
            "value": 0.9443,
            "uncertainty": 0.0072,
        },
        {
            "photon_eta": [1.5, 6.0],
            "photon_r9": [0.9, 999.0],
            "value": 0.9947,
            "uncertainty": 0.0051,
        },
    ],
}


def photon_preselection_sf(events, central_only):
    required_fields = [("Photon", "eta"), ("Photon", "r9")]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    if missing_fields:
        message = (
            "[photon_systematics : photon_preselection_sf] The events array is missing the following fields: %s which are needed as inputs."
            % (str(missing_fields))
        )
        logger.exception(message)
        raise ValueError(message)

    variations = systematic_from_bins(
        bins=PHOTON_PRESELECTION_SFs,
        variables={"photon_eta": abs(events.Photon.eta), "photon_r9": events.Photon.r9},
        central_only=central_only,
    )

    return variations


########################
### Electron veto SF ###
########################


photon_electron_veto_sf_bins = {
    "2016": PHOTON_ELECTRON_VETO_SF_2016,
    "2017": PHOTON_ELECTRON_VETO_SF_2017,
    "2018": PHOTON_ELECTRON_VETO_SF_2018,
}


def photon_electron_veto_sf(events, central_only, year):
    required_fields = [("Photon", "eta"), ("Photon", "r9")]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    if missing_fields:
        message = (
            "[photon_systematics : photon_preselection_sf] The events array is missing the following fields: %s which are needed as inputs."
            % (str(missing_fields))
        )
        logger.exception(message)
        raise ValueError(message)

    variations = systematic_from_bins(
        bins=photon_electron_veto_sf_bins[year],
        variables={"photon_eta": abs(events.Photon.eta), "photon_r9": events.Photon.r9},
        central_only=central_only,
    )

    return variations


##################
### Trigger SF ###
##################
# Note: since the trigger sf applies separate SF for the lead/sublead photons,
# it's easiest to just cast this as an EventWeightSystematic (rather than ObjectWeightSystematic as we would typically do)
# and just multiply the lead/sublead variations manually by hand


lead_trigger_sf_bins = {
    "2016": LEAD_TRIGGER_SF_2016,
    "2017": LEAD_TRIGGER_SF_2017,
    "2018": LEAD_TRIGGER_SF_2018,
}
sublead_trigger_sf_bins = {
    "2016": SUBLEAD_TRIGGER_SF_2016,
    "2017": SUBLEAD_TRIGGER_SF_2017,
    "2018": SUBLEAD_TRIGGER_SF_2018,
}


def trigger_sf(events, central_only, year):
    required_fields = [("Photon", "eta"), ("Photon", "r9"), ("Photon", "pt")]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    if missing_fields:
        message = (
            "[photon_systematics : photon_preselection_sf] The events array is missing the following fields: %s which are needed as inputs."
            % (str(missing_fields))
        )
        logger.exception(message)
        raise ValueError(message)

    variations_lead = systematic_from_bins(
        bins=lead_trigger_sf_bins[year],
        variables={
            "photon_r9": events.LeadPhoton.r9,
            "photon_eta": abs(events.LeadPhoton.eta),
            "photon_pt": events.LeadPhoton.pt,
        },
        central_only=central_only,
    )

    variations_sublead = systematic_from_bins(
        bins=sublead_trigger_sf_bins[year],
        variables={
            "photon_r9": events.SubleadPhoton.r9,
            "photon_eta": abs(events.SubleadPhoton.eta),
            "photon_pt": events.SubleadPhoton.pt,
        },
        central_only=central_only,
    )

    # Multiply up/down/central variations together, following this treatment in flashgg:
    # https://github.com/cms-analysis/flashgg/blob/1453740b1e4adc7184d5d8aa8a981bdb6b2e5f8e/Systematics/interface/DiPhotonFromSinglePhotonViewBase.h#L85-L87
    variations = {}
    for key in variations_lead.keys():
        variations[key] = variations_lead[key] * variations_sublead[key]

    return variations


############
### FNUF ###
############


fnuf_bins = {"2016": FNUF_2016, "2017": FNUF_2017, "2018": FNUF_2018}


def fnuf_unc(events, year, nominal_only, modify_nominal, loc="all"):
    """ """
    required_fields = [("Photon", "eta"), ("Photon", "r9")]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    if missing_fields:
        message = (
            "[photon_systematics : photon_preselection_sf] The events array is missing the following fields: %s which are needed as inputs."
            % (str(missing_fields))
        )
        logger.exception(message)
        raise ValueError(message)

    photons = events.Photon

    if loc == "all":
        mask = photons.pt > 0
    elif loc == "eb":
        mask = photons.isScEtaEB == True
    elif loc == "ee":
        mask = photons.isScEtaEE == True

    variations = ic_systematic_from_bins(
        bins=fnuf_bins[year],
        variables={"photon_eta": abs(events.Photon.eta), "photon_r9": events.Photon.r9},
        branch=photons.pt,
        nominal_only=nominal_only,
        modify_nominal=modify_nominal,
        mask=mask,
    )

    return variations


################
### Material ###
################


material_bins = {"2016": MATERIAL_2016, "2017": MATERIAL_2017, "2018": MATERIAL_2018}


def material_unc(events, year, nominal_only, modify_nominal, loc="all"):
    photons = events.Photon

    if loc == "all":
        mask = photons.pt > 0
    elif loc == "central_barrel":
        mask = abs(photons.eta) <= 1.0
    elif loc == "outer_barrel":
        mask = (abs(photons.eta) > 1.0) & (abs(photons.eta) <= 1.5)
    elif loc == "forward":
        mask = abs(photons.eta) > 1.5

    variations = ic_systematic_from_bins(
        bins=material_bins[year],
        variables={"photon_eta": abs(events.Photon.eta), "photon_r9": events.Photon.r9},
        branch=photons.pt,
        nominal_only=nominal_only,
        modify_nominal=modify_nominal,
        mask=mask,
    )

    return variations


def dummy_photon_pt_syst(events):
    photons = events.Photon

    variations = {}
    variations["up"] = photons.pt + awkward.ones_like(photons.pt)
    variations["down"] = photons.pt - awkward.ones_like(photons.pt)
    return variations


#######################
### Photon MC Smear ###
#######################


def photon_mc_smear(events, r9, loc):
    photons = events.Photon

    # Split into high/low r9 and EE/EB
    if r9 == "high":
        r9_cut = photons.r9 > 0.94
    elif r9 == "low":
        r9_cut = photons.r9 <= 0.94

    if loc == "eb":
        loc_cut = photons.isScEtaEB == True
    elif loc == "ee":
        loc_cut = photons.isScEtaEE == True

    mask = r9_cut & loc_cut

    variations = {}
    variations["up"] = awkward.where(mask, photons.pt + photons.dEsigmaUp, photons.pt)
    variations["down"] = awkward.where(
        mask, photons.pt + photons.dEsigmaDown, photons.pt
    )

    return variations
