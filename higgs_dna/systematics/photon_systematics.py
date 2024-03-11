import numpy as np
import awkward as ak
import correctionlib
import os
from copy import deepcopy
import logging

logger = logging.getLogger(__name__)


# first dummy, keeping it at this point as reference for even simpler implementations
def photon_pt_scale_dummy(pt, **kwargs):
    return (1.0 + np.array([0.05, -0.05], dtype=np.float32)) * pt[:, None]


# Not nice but working: if the functions are called in the base processor by Photon.add_systematic(... "what"="pt"...), the pt is passed to the function as first argument.
# I need the full events here, so I pass in addition the events. Seems to only work if it is explicitly a function of pt, but I might be missing something. Open for better solutions.
def Scale(pt, events, year="2022postEE", is_correction=True):
    """
    Applies the photon pt scale corrections (use on data!) and corresponding uncertainties (on MC!).
    JSONs need to be pulled first with scripts/pull_files.py
    """

    # for later unflattening:
    counts = ak.num(events.Photon.pt)

    run = ak.flatten(ak.broadcast_arrays(events.run, events.Photon.pt)[0])
    gain = ak.flatten(events.Photon.seedGain)
    eta = ak.flatten(events.Photon.ScEta)
    r9 = ak.flatten(events.Photon.r9)
    _pt = ak.flatten(events.Photon.pt)

    if year == "2016preVFP":
        path_json = os.path.join(os.path.dirname(__file__), 'JSONs/scaleAndSmearing/EGM_ScaleUnc_2016preVFP.json')
        evaluator = correctionlib.CorrectionSet.from_file(path_json)["UL-EGM_ScaleUnc"]
    elif year == "2016postVFP":
        path_json = os.path.join(os.path.dirname(__file__), 'JSONs/scaleAndSmearing/EGM_ScaleUnc_2016postVFP.json')
        evaluator = correctionlib.CorrectionSet.from_file(path_json)["UL-EGM_ScaleUnc"]
    elif year == "2017":
        path_json = os.path.join(os.path.dirname(__file__), 'JSONs/scaleAndSmearing/EGM_ScaleUnc_2017.json')
        evaluator = correctionlib.CorrectionSet.from_file(path_json)["UL-EGM_ScaleUnc"]
    elif year == "2018":
        path_json = os.path.join(os.path.dirname(__file__), 'JSONs/scaleAndSmearing/EGM_ScaleUnc_2018.json')
        evaluator = correctionlib.CorrectionSet.from_file(path_json)["UL-EGM_ScaleUnc"]
    elif year == "2022preEE":
        path_json = os.path.join(os.path.dirname(__file__), 'JSONs/scaleAndSmearing/SS_Rereco2022BCD.json')
        evaluator = correctionlib.CorrectionSet.from_file(path_json)["2022Re-recoBCD_ScaleJSON"]
    elif year == "2022postEE":
        path_json = os.path.join(os.path.dirname(__file__), 'JSONs/scaleAndSmearing/SS_RerecoE_PromptFG_2022.json')
        evaluator = correctionlib.CorrectionSet.from_file(path_json)["2022Re-recoE+PromptFG_ScaleJSON"]
    else:
        print("\n WARNING: there are only scale corrections for the year strings [\"2016preVFP\", \"2016postVFP\", \"2017\", \"2018\", \"2022preEE\", \"2022postEE\"]! \n Exiting. \n")
        exit()

    if is_correction:

        if year in ["2016preVFP", "2016postVFP", "2017", "2018"]:
            # the correction is already applied for Run 2
            logger.info("the scale correction for Run 2  MC is alsready applied in nAOD, nothing to be done")
        else:
            correction = evaluator.evaluate("total_correction", gain, run, eta, r9, _pt)
            pt_corr = _pt * correction

            corrected_photons = deepcopy(events.Photon)
            pt_corr = ak.unflatten(pt_corr, counts)
            corrected_photons["pt"] = pt_corr

            events.Photon = corrected_photons

        return events

    else:

        if year in ["2016preVFP", "2016postVFP", "2017", "2018"]:
            # the uncertainty is applied in reverse because the correction is meant for data as I understand fro EGM instructions here: https://cms-talk.web.cern.ch/t/pnoton-energy-corrections-in-nanoaod-v11/34327/2
            uncertainty_up = evaluator.evaluate(year, "scaledown", eta, gain)
            uncertainty_down = evaluator.evaluate(year, "scaleup", eta, gain)

            corr_up_variation = uncertainty_up
            corr_down_variation = uncertainty_down

        else:
            correction = evaluator.evaluate("total_correction", gain, run, eta, r9, _pt)
            uncertainty = evaluator.evaluate("total_uncertainty", gain, run, eta, r9, _pt)

            # divide by correction since it is already applied before
            corr_up_variation = (correction + uncertainty) / correction
            corr_down_variation = (correction - uncertainty) / correction

        # coffea does the unflattenning step itself and sets this value as pt of the up/down variations
        return np.concatenate((corr_up_variation.reshape(-1,1), corr_down_variation.reshape(-1,1)), axis=1) * _pt[:, None]


def Smearing(pt, events, year="2022postEE", is_correction=True):
    """
    Applies the photon smearing corrections and corresponding uncertainties (on MC!).
    JSON needs to be pulled first with scripts/pull_files.py
    """

    # for later unflattening:
    counts = ak.num(events.Photon.pt)

    eta = ak.flatten(events.Photon.ScEta)
    r9 = ak.flatten(events.Photon.r9)
    _pt = ak.flatten(events.Photon.pt)

    # we need reproducible random numbers since in the systematics call, the previous correction needs to be cancelled out
    rng = np.random.default_rng(seed=125)

    if year == "2022preEE":
        path_json = os.path.join(os.path.dirname(__file__), 'JSONs/scaleAndSmearing/SS_Rereco2022BCD.json')
        evaluator = correctionlib.CorrectionSet.from_file(path_json)["2022Re-recoBCD_SmearingJSON"]
    elif year == "2022postEE":
        path_json = os.path.join(os.path.dirname(__file__), 'JSONs/scaleAndSmearing/SS_RerecoE_PromptFG_2022.json')
        evaluator = correctionlib.CorrectionSet.from_file(path_json)["2022Re-recoE+PromptFG_SmearingJSON"]
    elif year in ["2016preVFP", "2016postVFP", "2017", "2018"]:
        logger.info("the systematic variations are taken directly from the dedicated nAOD branches Photon.dEsigmaUp and Photon.dEsigmaDown")
    else:
        print("\n WARNING: the correction for the selected year is not implemented yet! Valid year tags are [\"2016preVFP\", \"2016postVFP\", \"2017\", \"2018\", \"2022preEE\", \"2022postEE\"] \n Exiting. \n")
        exit()

    if is_correction:

        if year in ["2016", "2016preVFP", "2016postVFP", "2017", "2018"]:
            logger.info("the smearing correction for Run 2 MC is already applied in nAOD")
        else:
            # In theory, the energy should be smeared and not the pT, see: https://mattermost.web.cern.ch/cmseg/channels/egm-ss/6mmucnn8rjdgt8x9k5zaxbzqyh
            # However, there is a linear proportionality between pT and E: E = pT * cosh(eta)
            # Because of that, applying the correction to pT and E is equivalent (since eta does not change)
            # Energy is provided as a LorentzVector mixin, so we choose to correct pT
            # Also holds true for the scale part
            rho = evaluator.evaluate("rho", eta, r9)
            smearing = rng.normal(loc=1., scale=rho)
            pt_corr = _pt * smearing
            corrected_photons = deepcopy(events.Photon)
            pt_corr = ak.unflatten(pt_corr, counts)
            rho_corr = ak.unflatten(rho, counts)

            # If it is data, dont perform the pt smearing, only save the std of the gaussian for each event!
            try:
                events.GenIsolatedPhoton  # this operation is here because if there is no "events.GenIsolatedPhoton" field on data, an error will be thrown and we go to the except - so we dont smear the data pt spectrum
                corrected_photons["pt"] = pt_corr
            except:
                pass

            corrected_photons["rho_smear"] = rho_corr

            events.Photon = corrected_photons  # why does this work? Why do we not need events['Photon'] to assign?

        return events

    else:

        if year in ["2016", "2016preVFP", "2016postVFP", "2017", "2018"]:
            # the correction is already applied for Run 2
            dEsigmaUp = ak.flatten(events.Photon.dEsigmaUp)
            dEsigmaDown = ak.flatten(events.Photon.dEsigmaDown)
            logger.info(f"{dEsigmaUp}, {events.Photon.dEsigmaUp}")

            # the correction is given as additive factor for the Energy (Et_smear_up = Et + abs(dEsigmaUp)) so we have to convert it before applying it to the Pt
            # for EGM instruction on how to calculate uncertainties I link to this CMSTalk post https://cms-talk.web.cern.ch/t/pnoton-energy-corrections-in-nanoaod-v11/34327/2
            smearing_up = (_pt + abs(dEsigmaUp) / np.cosh(eta))
            smearing_down = (_pt - abs(dEsigmaDown) / np.cosh(eta))

            # divide by correction since it is already applied before we also divide for the Pt because it is later multipled when passing it to coffea, to be compatible with Run 3 calculation from json.
            # I convert it to numpy because ak.Arrays don't have the .reshape method needed further on.
            corr_up_variation = (smearing_up / _pt).to_numpy()
            corr_down_variation = (smearing_down / _pt).to_numpy()

        else:
            rho = evaluator.evaluate("rho", eta, r9)
            # produce the same numbers as in correction step
            smearing = rng.normal(loc=1., scale=rho)

            err_rho = evaluator.evaluate("err_rho", eta, r9)
            rho_up = rho + err_rho
            rho_down = rho - err_rho
            smearing_up = rng.normal(loc=1., scale=rho_up)
            smearing_down = rng.normal(loc=1., scale=rho_down)

            # divide by correction since it is already applied before
            corr_up_variation = smearing_up / smearing
            corr_down_variation = smearing_down / smearing

        # coffea does the unflattenning step itself and sets this value as pt of the up/down variations
        return np.concatenate((corr_up_variation.reshape(-1,1), corr_down_variation.reshape(-1,1)), axis=1) * _pt[:, None]


def energyErrShift(energyErr, events, year="2022postEE", is_correction=True):
    # See also https://indico.cern.ch/event/1131803/contributions/4758593/attachments/2398621/4111806/Hgg_Differentials_Approval_080322.pdf#page=47
    if is_correction:
        return events
    else:
        _energyErr = ak.flatten(events.Photon.energyErr)
        uncertainty_up = np.ones(len(_energyErr)) * 1.05
        uncertainty_dn = np.ones(len(_energyErr)) * 0.95
        return (
            np.concatenate(
                (uncertainty_up.reshape(-1, 1), uncertainty_dn.reshape(-1, 1)), axis=1
            )
            * _energyErr[:, None]
        )


# Not nice but working: if the functions are called in the base processor by Photon.add_systematic(... "what"="pt"...), the pt is passed to the function as first argument.
# I need the full events here, so I pass in addition the events. Seems to only work if it is explicitly a function of pt, but I might be missing something. Open for better solutions.
def FNUF(pt, events, year="2017", is_correction=True):
    """
    ---This is an implementation of the FNUF uncertainty copied from flashgg,
    --- Preliminary JSON (run2 I don't know if this needs to be changed) file created with correctionlib starting from flashgg: https://github.com/cms-analysis/flashgg/blob/2677dfea2f0f40980993cade55144636656a8a4f/Systematics/python/flashggDiPhotonSystematics2017_Legacy_cfi.py
    Applies the photon pt and energy scale corrections and corresponding uncertainties.
    To be checked by experts
    """

    # for later unflattening:
    counts = ak.num(events.Photon.pt)
    eta = ak.flatten(events.Photon.ScEta)
    r9 = ak.flatten(events.Photon.r9)
    _energy = ak.flatten(events.Photon.energy)
    _pt = ak.flatten(events.Photon.pt)

    # era/year defined as parameter of the function, only 2017 is implemented up to now
    avail_years = ["2016", "2016preVFP", "2016postVFP", "2017", "2018"]
    if year not in avail_years:
        print(f"\n WARNING: only FNUF corrections for the year strings {avail_years} are already implemented! \n Exiting. \n")
        exit()
    elif "2016" in year:
        year = "2016"

    jsonpog_file = os.path.join(os.path.dirname(__file__), f"JSONs/FNUF/{year}/FNUF_{year}.json")
    evaluator = correctionlib.CorrectionSet.from_file(jsonpog_file)["FNUF"]

    if is_correction:
        correction = evaluator.evaluate("nominal", eta, r9)
        corr_energy = _energy * correction
        corr_pt = _pt * correction

        corrected_photons = deepcopy(events.Photon)
        corr_energy = ak.unflatten(corr_energy, counts)
        corr_pt = ak.unflatten(corr_pt, counts)
        corrected_photons["energy"] = corr_energy
        corrected_photons["pt"] = corr_pt
        events.Photon = corrected_photons

        return events

    else:
        correction = evaluator.evaluate("nominal", eta, r9)
        # When creating the JSON I already included added the variation to the returned value
        uncertainty_up = evaluator.evaluate("up", eta, r9) / correction
        uncertainty_dn = evaluator.evaluate("down", eta, r9) / correction
        # coffea does the unflattenning step itself and sets this value as pt of the up/down variations
        return (
            np.concatenate(
                (uncertainty_up.reshape(-1, 1), uncertainty_dn.reshape(-1, 1)), axis=1
            )
            * _pt[:, None]
        )


# Same old same old, just reiterated: if the functions are called in the base processor by Photon.add_systematic(... "what"="pt"...), the pt is passed to the function as first argument.
# Open for better solutions.
def ShowerShape(pt, events, year="2017", is_correction=True):
    """
    ---This is an implementation of the ShowerShape uncertainty copied from flashgg,
    --- Preliminary JSON (run2 I don't know if this needs to be changed) file created with correctionlib starting from flashgg: https://github.com/cms-analysis/flashgg/blob/2677dfea2f0f40980993cade55144636656a8a4f/Systematics/python/flashggDiPhotonSystematics2017_Legacy_cfi.py
    Applies the photon pt and energy scale corrections and corresponding uncertainties (only on the pt because it is what is used in selection).
    To be checked by experts
    """

    # for later unflattening:
    counts = ak.num(events.Photon.pt)
    eta = ak.flatten(events.Photon.ScEta)
    r9 = ak.flatten(events.Photon.r9)
    _energy = ak.flatten(events.Photon.energy)
    _pt = ak.flatten(events.Photon.pt)

    # era/year defined as parameter of the function
    avail_years = ["2016", "2016preVFP", "2016postVFP", "2017", "2018"]
    if year not in avail_years:
        print(f"\n WARNING: only ShowerShape corrections for the year strings {avail_years} are already implemented! \n Exiting. \n")
        exit()
    elif "2016" in year:
        year = "2016"

    jsonpog_file = os.path.join(os.path.dirname(__file__), f"JSONs/ShowerShape/{year}/ShowerShape_{year}.json")
    evaluator = correctionlib.CorrectionSet.from_file(jsonpog_file)["ShowerShape"]

    if is_correction:
        correction = evaluator.evaluate("nominal", eta, r9)
        corr_energy = _energy * correction
        corr_pt = _pt * correction

        corrected_photons = deepcopy(events.Photon)
        corr_energy = ak.unflatten(corr_energy, counts)
        corr_pt = ak.unflatten(corr_pt, counts)
        corrected_photons["energy"] = corr_energy
        corrected_photons["pt"] = corr_pt
        events.Photon = corrected_photons

        return events

    else:
        correction = evaluator.evaluate("nominal", eta, r9)
        # When creating the JSON I already included added the variation to the returned value
        uncertainty_up = evaluator.evaluate("up", eta, r9) / correction
        uncertainty_dn = evaluator.evaluate("down", eta, r9) / correction
        # coffea does the unflattenning step itself and sets this value as pt of the up/down variations
        return (
            np.concatenate(
                (uncertainty_up.reshape(-1, 1), uncertainty_dn.reshape(-1, 1)), axis=1
            )
            * _pt[:, None]
        )


def Material(pt, events, year="2017", is_correction=True):
    """
    ---This is an implementation of the Material uncertainty copied from flashgg,
    --- JSON file for run2 created with correctionlib starting from flashgg: https://github.com/cms-analysis/flashgg/blob/2677dfea2f0f40980993cade55144636656a8a4f/Systematics/python/flashggDiPhotonSystematics2017_Legacy_cfi.py
    Applies the photon pt and energy scale corrections and corresponding uncertainties.
    To be checked by experts
    """

    # for later unflattening:
    counts = ak.num(events.Photon.pt)
    eta = ak.flatten(abs(events.Photon.ScEta))
    r9 = ak.flatten(events.Photon.r9)
    _energy = ak.flatten(events.Photon.energy)
    _pt = ak.flatten(events.Photon.pt)

    # era/year defined as parameter of the function, only 2017 is implemented up to now
    avail_years = ["2016", "2016preVFP", "2016postVFP", "2017", "2018", "2022preEE", "2022postEE"]
    if year not in avail_years:
        print(f"\n WARNING: only eVetoSF corrections for the year strings {avail_years} are already implemented! \n Exiting. \n")
        exit()
    elif "2016" in year:
        year = "2016"
    # use Run 2 files also for Run 3, preliminary
    elif year in ["2022preEE", "2022postEE"]:
        year = "2017"

    jsonpog_file = os.path.join(os.path.dirname(__file__), f"JSONs/Material/{year}/Material_{year}.json")
    evaluator = correctionlib.CorrectionSet.from_file(jsonpog_file)["Material"]

    if is_correction:
        correction = evaluator.evaluate("nominal", eta, r9)
        corr_energy = _energy * correction
        corr_pt = _pt * correction

        corrected_photons = deepcopy(events.Photon)
        corr_energy = ak.unflatten(corr_energy, counts)
        corr_pt = ak.unflatten(corr_pt, counts)
        corrected_photons["energy"] = corr_energy
        corrected_photons["pt"] = corr_pt
        events.Photon = corrected_photons

        return events

    else:
        correction = evaluator.evaluate("nominal", eta, r9)
        # When creating the JSON I already included added the variation to the returned value
        uncertainty_up = evaluator.evaluate("up", eta, r9) / correction
        uncertainty_dn = evaluator.evaluate("down", eta, r9) / correction
        # coffea does the unflattenning step itself and sets this value as pt of the up/down variations
        return (
            np.concatenate(
                (uncertainty_up.reshape(-1, 1), uncertainty_dn.reshape(-1, 1)), axis=1
            )
            * _pt[:, None]
        )
