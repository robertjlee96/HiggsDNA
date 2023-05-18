import numpy as np
import awkward as ak
import os
import uproot


def add_pileup_weight(events: ak.Array, pileup_profile: str = None) -> ak.Array:
    """
    Parameter pileup_profile not yet used here. Idea: integrate in metaconditions json which file to read in this function.
    Reading privately produced file here.
    """
    path_pileup_profile = os.path.join(
        os.path.dirname(__file__),
        "../metaconditions/pileup_histograms/MyDataPileupHistogram_2022Full.root",
    )
    pileup_profile = uproot.open(path_pileup_profile)["pileup"]
    pileup_profile = pileup_profile.to_numpy()[0]
    pileup_profile /= pileup_profile.sum()

    pileup_MC = np.histogram(ak.to_numpy(events.Pileup.nPU), bins=100, range=(0, 100))[0].astype("float64")
    # avoid division by zero later
    pileup_MC[pileup_MC == 0.] = 1
    pileup_MC /= pileup_MC.sum()

    pileup_correction = pileup_profile / pileup_MC
    # remove large MC reweighting factors to prevent artifacts
    pileup_correction[pileup_correction > 10] = 1

    weight_pileup = pileup_correction[ak.to_numpy(events.Pileup.nPU)]
    events["weight_pileup"] = weight_pileup

    return events
