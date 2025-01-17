from .photon_systematics import (
    photon_pt_scale_dummy,
    Scale,
    Smearing,
    energyErrShift,
    FNUF,
    ShowerShape,
    Material,
)
from .event_weight_systematics import (
    Pileup,
    SF_photon_ID,
    LooseMvaSF,
    ElectronVetoSF,
    PreselSF,
    TriggerSF,
    NNLOPS,
    AlphaS,
    PartonShower,
    cTagSF,
    Zpt,
)
from .jet_systematics import (
    jet_pt_scale_dummy,
    JERC_jet,
)
from .jet_systematics_json import jerc_jet
from functools import partial
import logging

logger = logging.getLogger(__name__)

# using add_systematic function of coffea.nanoevents.methods.nanoaod objects as Photon to store object systematics in addition to nominal objects
object_systematics = {
    "PhotonPtScale_dummy": {
        "object": "Photon",
        "args": {
            "kind": "UpDownSystematic",
            "what": "pt",
            "varying_function": photon_pt_scale_dummy,
        },
    },
    "Scale": {
        "object": "Photon",
        "args": {
            "kind": "UpDownSystematic",
            "what": "pt",
            "varying_function": partial(Scale, is_correction=False),
        },
    },
    "Smearing": {
        "object": "Photon",
        "args": {
            "kind": "UpDownSystematic",
            "what": "pt",
            "varying_function": partial(Smearing, is_correction=False),
        },
    },
    "energyErrShift": {
        "object": "Photon",
        "args": {
            "kind": "UpDownSystematic",
            "what": "energyErr",
            "varying_function": partial(energyErrShift, is_correction=False),
        },
    },
    "FNUF": {
        "object": "Photon",
        "args": {
            "kind": "UpDownSystematic",
            "what": "pt",
            "varying_function": partial(FNUF, year="2017", is_correction=False),
        },
    },
    "ShowerShape": {
        "object": "Photon",
        "args": {
            "kind": "UpDownSystematic",
            "what": "pt",
            "varying_function": partial(ShowerShape, year="2017", is_correction=False),
        },
    },
    "Material": {
        "object": "Photon",
        "args": {
            "kind": "UpDownSystematic",
            "what": "pt",
            "varying_function": partial(Material, year="2017", is_correction=False),
        },
    },
    "JetPtScale_dummy": {
        "object": "Jet",
        "args": {
            "kind": "UpDownSystematic",
            "what": "pt",
            "varying_function": partial(
                jet_pt_scale_dummy, year=None, is_correction=False
            ),
        },
    },
    "JES": {
        "object": "Jet",
        "args": {
            "kind": "UpDownSystematic",
            "what": "pt",
            "varying_function": partial(
                JERC_jet,
                year="2022postEE",
                skip_JER=False,
                skip_JEC=False,
                is_correction=False,
            ),
        },
    },
}

# functions correcting nominal object quantities to be placed here
# dict containing "name": varying_function
object_corrections = {
    "Scale": partial(Scale, pt=None, is_correction=True),
    "Smearing": partial(Smearing, pt=None, is_correction=True),
    "energyErrShift": partial(energyErrShift, energyErr=None, is_correction=True),
    "FNUF": partial(FNUF, pt=None, is_correction=True),
    "ShowerShape": partial(ShowerShape, pt=None, is_correction=True),
    "Material": partial(Material, pt=None, is_correction=True),
    "JetPtScale_dummy": partial(
        jet_pt_scale_dummy, pt=None, year=None, is_correction=True
    ),
    "JES": partial(
        JERC_jet,
        pt=None,
        year="2022postEE",
        skip_JER=False,
        skip_JEC=False,
        is_correction=True,
    ),
}

# functions adding systematic variations to event weights to be placed here
# dict containing "name": varying_function
weight_systematics = {
    "Pileup": partial(Pileup, is_correction=False),
    "SF_photon_ID": partial(SF_photon_ID, is_correction=False),
    "LooseMvaSF": partial(LooseMvaSF, is_correction=False),
    "ElectronVetoSF": partial(ElectronVetoSF, is_correction=False),
    "PreselSF": partial(PreselSF, is_correction=False),
    "TriggerSF": partial(TriggerSF, is_correction=False),
    "cTagSF": partial(cTagSF, is_correction=False),
    "AlphaS": partial(AlphaS),
    "PartonShower": partial(PartonShower),
    "LHEScale": None,
    "LHEPdf": None,
    "Zpt": partial(Zpt),
}

# functions correcting nominal event weights to be placed here
# dict containing "name": varying_function
weight_corrections = {
    "Pileup": partial(Pileup, is_correction=True),
    "SF_photon_ID": partial(SF_photon_ID, is_correction=True),
    "LooseMvaSF": partial(LooseMvaSF, is_correction=True),
    "ElectronVetoSF": partial(ElectronVetoSF, is_correction=True),
    "PreselSF": partial(PreselSF, is_correction=True),
    "TriggerSF": partial(TriggerSF, is_correction=True),
    "cTagSF": partial(cTagSF, is_correction=True),
    "NNLOPS": partial(NNLOPS, is_correction=True),
    "Zpt": partial(Zpt, is_correction=True),
}


def check_corr_syst_combinations(corrections_dict, systematics_dict, logger):
    """
    This function is a sanity check for the choice of systematics and corrections which the user wants to process.
    It ensures that systematic variations of a correction can only be processed when the correction itself is applied.
    """
    for dataset in systematics_dict.keys():
        for chosen_syst in systematics_dict[dataset]:
            if (
                chosen_syst in weight_corrections.keys()
                and chosen_syst not in corrections_dict[dataset]
            ) or (
                chosen_syst in object_corrections.keys()
                and chosen_syst not in corrections_dict[dataset]
            ):
                # scale unc. will be applied to MC while the correction is applied to data. Exception.
                if (
                    "scale" in chosen_syst.lower()
                    or "jec" in chosen_syst.lower()
                    or "jer" in chosen_syst.lower()
                ):
                    continue
                logger.info(
                    f"Requested to evaluate systematic variation {chosen_syst} for dataset {dataset} without applying the corresponding correction. \nThis is not intended.\nExiting."
                )
                exit()


def add_jme_corr_syst(corrections_dict, systematics_dict, logger):
    corrections_dict.update(
        {
            # jerc for MC
            "jec_jet": partial(jerc_jet, pt=None, apply_jec=True),
            "jec_jet_syst": partial(jerc_jet, pt=None, apply_jec=True, jec_syst=True),
            "jec_jet_regrouped_syst": partial(
                jerc_jet, pt=None, apply_jec=True, jec_syst=True, split_jec_syst=True
            ),
            "jerc_jet": partial(jerc_jet, pt=None, apply_jec=True, apply_jer=True),
            "jerc_jet_syst": partial(
                jerc_jet,
                pt=None,
                apply_jec=True,
                jec_syst=True,
                apply_jer=True,
                jer_syst=True,
            ),
            "jerc_jet_regrouped_syst": partial(
                jerc_jet,
                pt=None,
                apply_jec=True,
                jec_syst=True,
                split_jec_syst=True,
                apply_jer=True,
                jer_syst=True,
            ),
            # jec for data: Usually corrections on Data innecessary
            # Use jec corrections with Era to re-do jec for data
            "jec_RunA": partial(jerc_jet, pt=None, era="RunA", level="L1L2L3Res"),
            "jec_RunB": partial(jerc_jet, pt=None, era="RunB", level="L1L2L3Res"),
            "jec_RunC": partial(jerc_jet, pt=None, era="RunC", level="L1L2L3Res"),
            "jec_RunD": partial(jerc_jet, pt=None, era="RunD", level="L1L2L3Res"),
            "jec_RunE": partial(jerc_jet, pt=None, era="RunE", level="L1L2L3Res"),
            "jec_RunF": partial(jerc_jet, pt=None, era="RunF", level="L1L2L3Res"),
            "jec_RunG": partial(jerc_jet, pt=None, era="RunG", level="L1L2L3Res"),
            "jec_RunH": partial(jerc_jet, pt=None, era="RunH", level="L1L2L3Res"),
        }
    )
    logger.info(
        f"""_summary_

    Available correction keys:
    {corrections_dict.keys()}
    Available systematic keys:
    {systematics_dict.keys()}
    """
    )
    return corrections_dict, systematics_dict


object_corrections, object_systematics = add_jme_corr_syst(
    object_corrections, object_systematics, logger
)
