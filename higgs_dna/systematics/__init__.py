from .photon_systematics import photon_pt_scale_dummy, Scale, FNUF, ShowerShape
from .event_weight_systematics import SF_photon_ID, LooseMvaSF, NNLOPS
from functools import partial

# using add_systematic function of coffea.nanoevents.methods.nanoaod objects as Photon to store object systematics in addition to nominal objects
object_systematics = {
    "PhotonPtScale_dummy": {
        "object": "Photon",
        "args": {
            "kind": "UpDownSystematic",
            "what": "pt",
            "varying_function": photon_pt_scale_dummy,
        }
    },
    "Scale_2016postVFP": {
        "object": "Photon",
        "args": {
            "kind": "UpDownSystematic",
            "what": "pt",
            "varying_function": partial(Scale, year="2016postVFP", is_correction=False),
        }
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
}

# functions correcting nominal object quantities to be placed here
# dict containing "name": varying_function
object_corrections = {
    "Scale_2016postVFP": partial(Scale, pt=None,year="2016postVFP", is_correction=True),
    "FNUF": partial(FNUF, pt=None, year="2017", is_correction=True),
    "ShowerShape": partial(ShowerShape, pt=None, year="2017", is_correction=True),
}

# functions adding systematic variations to event weights to be placed here
# dict containing "name": varying_function
weight_systematics = {
    "SF_photon_ID": partial(SF_photon_ID, is_correction=False),
    "LooseMvaSF": partial(LooseMvaSF, is_correction=False),
}

# functions correcting nominal event weights to be placed here
# dict containing "name": varying_function
weight_corrections = {
    "SF_photon_ID": partial(SF_photon_ID, is_correction=True),
    "LooseMvaSF": partial(LooseMvaSF, is_correction=True),
    "NNLOPS": partial(NNLOPS, is_correction=True),
}


def check_corr_syst_combinations(corrections_dict, systematics_dict, logger):
    """
    This function is a sanity check for the choice of systematics and corrections which the user wants to process.
    It ensures that systematic variations of a correction can only be processed when the correction itself is applied.
    """
    for dataset in systematics_dict.keys():
        for chosen_syst in systematics_dict[dataset]:
            if (chosen_syst in weight_corrections.keys() and chosen_syst not in corrections_dict[dataset]) or (chosen_syst in object_corrections.keys() and chosen_syst not in corrections_dict[dataset]):
                logger.info(
                    f"Requested to evaluate systematic variation {chosen_syst} for dataset {dataset} without applying the corresponding correction. \nThis is not intended.\nExiting."
                )
                exit()
