import numpy as np
import json
import os
from scipy.interpolate import interp1d
import correctionlib


def SF_photon_ID(photons, weights, year="2017", WP="Loose", is_correction=True, **kwargs):
    """
    ---This is a dummy, meant to be replaced by the Run-3 photon ID continuous SF later, but shows how it can be implemented.---
    Applies the photon ID scale-factor and corresponding uncertainties
    as specified in https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration.
    **kwargs used here since I pass the full events object to these functions as some other systematics might need more information than only the Photon container.
    """

    # have to think about how to specify era/year, using 2017 test-wise here, defined as parameter of the function
    jsonpog_file = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2017_UL/photon.json.gz"
    evaluator = correctionlib.CorrectionSet.from_file(jsonpog_file)["UL-Photon-ID-SF"]

    if is_correction:
        # only calculate correction to nominal weight
        sf_lead = evaluator.evaluate(
            year, "sf", "Loose", photons["pho_lead"].eta, photons["pho_lead"].pt
        )
        sf_sublead = evaluator.evaluate(
            year, "sf", "Loose", photons["pho_sublead"].eta, photons["pho_sublead"].pt
        )
        sf = sf_lead * sf_sublead

        sfup, sfdown = None, None

    else:
        # only calculate systs
        sf = np.ones(len(weights._weight))
        sf_lead = evaluator.evaluate(
            year, "sf", "Loose", photons["pho_lead"].eta, photons["pho_lead"].pt
        )
        sf_sublead = evaluator.evaluate(
            year, "sf", "Loose", photons["pho_sublead"].eta, photons["pho_sublead"].pt
        )
        _sf = sf_lead * sf_sublead

        sfup_lead = evaluator.evaluate(
            year, "sfup", "Loose", photons["pho_lead"].eta, photons["pho_lead"].pt
        )
        sfup_sublead = evaluator.evaluate(
            year, "sfup", "Loose", photons["pho_sublead"].eta, photons["pho_sublead"].pt
        )
        # Temporary fix: the .weight function of coffea.processor.Weights multiplies the systematic variations by the nominal SF itself
        # so divide by it to remove this additional factor
        sfup = sfup_lead * sfup_sublead / _sf

        sfdown_lead = evaluator.evaluate(
            year, "sfdown", "Loose", photons["pho_lead"].eta, photons["pho_lead"].pt
        )
        sfdown_sublead = evaluator.evaluate(
            year, "sfdown", "Loose", photons["pho_sublead"].eta, photons["pho_sublead"].pt,
        )
        sfdown = sfdown_lead * sfdown_sublead / _sf

    weights.add(name="SF_photon_ID", weight=sf, weightUp=sfup, weightDown=sfdown)

    return weights


def LooseMvaSF(photons, weights, year="2017", WP="Loose", is_correction=True, **kwargs):
    """
    LooseMvaSF: correction to the event weight on a per photon level, impacting one of the high importance input variable of the DiphotonBDT, binned in eta and r9.
    for original implementation look at: https://github.com/cms-analysis/flashgg/blob/2677dfea2f0f40980993cade55144636656a8a4f/Systematics/python/flashggDiPhotonSystematics2017_Legacy_cfi.py
    And for presentation on the study: https://indico.cern.ch/event/963617/contributions/4103623/attachments/2141570/3608645/Zee_Validation_UL2017_Update_09112020_Prasant.pdf
    """

    # have to think about how to specify era/year, using 2017 test-wise here, defined as parameter of the function
    jsonpog_file = os.path.join(os.path.dirname(__file__), "JSONs/LooseMvaSF.json")
    evaluator = correctionlib.CorrectionSet.from_file(jsonpog_file)["LooseMvaSF"]

    if is_correction:
        # only calculate correction to nominal weight
        sf_lead = evaluator.evaluate(
            "nominal", photons["pho_lead"].eta, photons["pho_lead"].r9
        )
        sf_sublead = evaluator.evaluate(
            "nominal", photons["pho_sublead"].eta, photons["pho_sublead"].r9
        )
        sf = sf_lead * sf_sublead

        sfup, sfdown = None, None

    else:
        # only calculate systs
        sf = np.ones(len(weights._weight))
        sf_lead = evaluator.evaluate(
            "nominal", photons["pho_lead"].eta, photons["pho_lead"].r9
        )
        sf_sublead = evaluator.evaluate(
            "nominal", photons["pho_sublead"].eta, photons["pho_sublead"].r9
        )
        _sf = sf_lead * sf_sublead

        sfup_lead = evaluator.evaluate(
            "up", photons["pho_lead"].eta, photons["pho_lead"].r9
        )
        sfup_sublead = evaluator.evaluate(
            "up", photons["pho_sublead"].eta, photons["pho_sublead"].r9
        )
        sfup = sfup_lead * sfup_sublead / _sf

        sfdown_lead = evaluator.evaluate(
            "down", photons["pho_lead"].eta, photons["pho_lead"].r9
        )
        sfdown_sublead = evaluator.evaluate(
            "down", photons["pho_sublead"].eta, photons["pho_sublead"].r9
        )
        sfdown = sfdown_lead * sfdown_sublead / _sf

    weights.add(name="LooseMvaSF", weight=sf, weightUp=sfup, weightDown=sfdown)

    return weights


def NNLOPS(events, dataset_name, weights, is_correction=True, generator="mcatnlo", **kwargs):
    """
    --- NNLOPS reweighting for ggH events to be applied to NLO Madgraph (and Powheg).
    Swap generator argument to 'powheg' if to be applied to powheg events
    Reweight event based on truth Higgs pt and number of jets, extracted from HTXS object
    Constructs njet-dependent linear splines based on input data, functions of Higgs pt
    Reweighting only applied if ggH is in the dataset name, otherwise get a scale factor of 1
    """
    json_file = os.path.join(os.path.dirname(__file__), "JSONs/NNLOPS_reweight.json")

    if is_correction:

        if "ggH" in dataset_name:

            # Extract NNLOPS weights from json file
            with open(json_file, "r") as jf:
                nnlops_reweight = json.load(jf)

            # Load reweight factors for specific generator
            nnlops_reweight = nnlops_reweight[generator]

            # Build linear splines for different njet bins
            spline_0jet = interp1d(nnlops_reweight['0jet']['pt'], nnlops_reweight['0jet']['weight'])
            spline_1jet = interp1d(nnlops_reweight['1jet']['pt'], nnlops_reweight['1jet']['weight'])
            spline_2jet = interp1d(nnlops_reweight['2jet']['pt'], nnlops_reweight['2jet']['weight'])
            spline_ge3jet = interp1d(nnlops_reweight['3jet']['pt'], nnlops_reweight['3jet']['weight'])

            # Load truth Higgs pt and njets (pt>30) from events
            higgs_pt = events.HTXS.Higgs_pt
            njets30 = events.HTXS.njets30

            # Extract scale factors from splines and mask for different jet bins
            # Define maximum pt values as interpolated splines only go up so far
            sf = (njets30 == 0) * spline_0jet(np.minimum(np.array(higgs_pt), 125.)) + \
                 (njets30 == 1) * spline_1jet(np.minimum(np.array(higgs_pt), 625.)) + \
                 (njets30 == 2) * spline_2jet(np.minimum(np.array(higgs_pt), 800.)) + \
                 (njets30 >= 3) * spline_ge3jet(np.minimum(np.array(higgs_pt), 925.))

        else:
            sf = np.ones(len(weights._weight))

    else:
        raise RuntimeError("NNLOPS reweighting is only a flat correction, not a systematic")

    weights.add("NNLOPS", sf, None, None)

    return weights
