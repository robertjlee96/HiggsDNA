import numpy as np
import awkward as ak
import correctionlib
import os
from copy import deepcopy


# first dummy, keeping it at this point as reference for even simpler implementations
def photon_pt_scale_dummy(pt, **kwargs):
    return (1.0 + np.array([0.05, -0.05], dtype=np.float32)) * pt[:, None]


# Not nice but working: if the functions are called in the base processor by Photon.add_systematic(... "what"="pt"...), the pt is passed to the function as first argument.
# I need the full events here, so I pass in addition the events. Seems to only work if it is explicitly a function of pt, but I might be missing something. Open for better solutions.
def Scale(pt, events, year="2016postVFP", is_correction=True):
    """
    ---This is a dummy, meant to be replaced by the Run-3 photon scale uncertainties later, which only shows how a scale uncertainty/correction can be implemented---
    --- Preliminary JSON file taken from https://github.com/a-kapoor/ScaleFactorsJSON/tree/master/2016postVFP_UL ---
    Applies the photon pt scale corrections and corresponding uncertainties.
    The official JSON to use should apper in https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration later.
    """

    # for later unflattening:
    counts = ak.num(events.Photon.pt)
    # dummy value to test on MC/data which are not 2016postVFP:
    run = 278769.
    # otherwise uncomment the following line:
    # run = ak.flatten(events.run)
    gain = ak.flatten(events.Photon.seedGain)
    eta = ak.flatten(events.Photon.eta)
    # this JSON does not have all r9 bins. Setting it to 0.8 here for tests
    r9 = 0.8
    # r9 = ak.flatten(events.Photon.r9)
    _pt = ak.flatten(events.Photon.pt)

    jsonpog_file = os.path.join(os.path.dirname(__file__), 'JSONs/SS_2016postVFP_preliminary.json')
    evaluator = correctionlib.CorrectionSet.from_file(jsonpog_file)["2016postVFP_ScaleJSON"]

    if is_correction:

        correction = evaluator.evaluate("total_correction", gain, run, eta, r9, _pt)
        pt_corr = _pt * correction

        corrected_photons = deepcopy(events.Photon)
        pt_corr = ak.unflatten(pt_corr, counts)
        corrected_photons["pt"] = pt_corr

        events.Photon = corrected_photons

        return events

    else:

        correction = evaluator.evaluate("total_correction", gain, run, eta, r9, _pt)
        uncertainty = evaluator.evaluate("total_uncertainty", gain, run, eta, r9, _pt)

        # divide by correction since it is already applied before
        corr_up_variation = (correction + uncertainty) / correction
        corr_down_variation = (correction - uncertainty) / correction

        # coffea does the unflattenning step itself and sets this value as pt of the up/down variations
        return np.concatenate((corr_up_variation.reshape(-1,1), corr_down_variation.reshape(-1,1)), axis=1) * _pt[:, None]
