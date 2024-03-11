import numpy as np
import json
import os
from scipy.interpolate import interp1d
import correctionlib
import awkward as ak
from higgs_dna.utils.misc_utils import choose_jet
import logging

logger = logging.getLogger(__name__)


def SF_photon_ID(
    photons, weights, year="2017", WP="Loose", is_correction=True, **kwargs
):
    """
    Applies the photon ID scale-factor and corresponding uncertainties for the customised cut on the EGamma MVA ID (Run 3)
    JLS removed the support for the EGamma MVA ID SFs from https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration for Run 2 for now as this is not commonly used in the Hgg group
    Take action yourself or contact us if you need those!
    """
    # era/year defined as parameter of the function
    avail_years = ["2022preEE", "2022postEE"]
    if year not in avail_years:
        print(f"\n WARNING: only scale corrections for the year strings {avail_years} are already implemented! \n Exiting. \n")
        print("If you need the SFs for the central Egamma MVA ID for Run 2 UL, take action yourself or contact us!")
        exit()

    if year == "2022preEE":
        json_file = os.path.join(os.path.dirname(__file__), "JSONs/SF_photon_ID/2022/PhotonIDMVA_2022PreEE.json")
    elif year == "2022postEE":
        json_file = os.path.join(os.path.dirname(__file__), "JSONs/SF_photon_ID/2022/PhotonIDMVA_2022PostEE.json")

    evaluator = correctionlib.CorrectionSet.from_file(json_file)["PhotonIDMVA_SF"]

    # In principle, we should use the fully correct formula https://indico.cern.ch/event/1360948/contributions/5783762/attachments/2788516/4870824/24_02_02_HIG-23-014_PreAppPres.pdf#page=7
    # However, if the SF is pt-binned, the approximation of the multiplication of the two SFs is fully exact
    if "2022" in year:
        if is_correction:
            # only calculate correction to nominal weight
            sf_lead = evaluator.evaluate(
                abs(photons["pho_lead"].ScEta), photons["pho_lead"].pt, "nominal"
            )
            sf_sublead = evaluator.evaluate(
                abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].pt, "nominal"
            )
            sf = sf_lead * sf_sublead

            sfup, sfdown = None, None

        else:
            # only calculate systs

            sf = np.ones(len(weights._weight))
            sf_lead = evaluator.evaluate(
                abs(photons["pho_lead"].ScEta), photons["pho_lead"].pt, "nominal"
            )
            sf_sublead = evaluator.evaluate(
                abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].pt, "nominal"
            )
            _sf = sf_lead * sf_sublead

            sf_unc_lead = evaluator.evaluate(
                abs(photons["pho_lead"].ScEta), photons["pho_lead"].pt, "uncertainty"
            )
            sf_unc_sublead = evaluator.evaluate(
                abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].pt, "uncertainty"
            )

            sfup = (sf_lead + sf_unc_lead) * (sf_sublead + sf_unc_sublead) / _sf

            sfdown = (sf_lead - sf_unc_lead) * (sf_sublead - sf_unc_sublead) / _sf

    weights.add(name="SF_photon_ID", weight=sf, weightUp=sfup, weightDown=sfdown)

    return weights


def Pileup(events, weights, year, is_correction=True, **kwargs):
    """
    Function to apply either the pileup correction to MC to make it match the pileup profile of a certain year/period,
    or the respective uncertainties.
    The parameter `year` needs to be specified as one of ["2022preEE", "2022postEE", "23preBPix", "23postBPix"] for Run-3 or ["2016preVFP", "2016postVFP", "2017", "2018"] for Run-2.
    By now, the Run-2 and Run-3 up to 2023D files are available from LUM POG in the correctionlib format...
    The pileup histos for Run-3 were produced by Junquan, the JSONs for Run-2 and Run-3 first need to be pulled with `scripts/pull_files.py`!
    """
    path_to_json = os.path.join(os.path.dirname(__file__), "../systematics/JSONs/pileup/pileup_{}.json.gz".format(year))
    if "16" in year:
        name = "Collisions16_UltraLegacy_goldenJSON"
    elif "17" in year:
        name = "Collisions17_UltraLegacy_goldenJSON"
    elif "18" in year:
        name = "Collisions18_UltraLegacy_goldenJSON"
    elif "22preEE" in year:
        name = "Collisions2022_355100_357900_eraBCD_GoldenJson"
    elif "22postEE" in year:
        name = "Collisions2022_359022_362760_eraEFG_GoldenJson"
    elif "23preBPix" in year:
        name = "Collisions2023_366403_369802_eraBC_GoldenJson"
    elif "23postBPix" in year:
        name = "Collisions2023_369803_370790_eraD_GoldenJson"

    evaluator = correctionlib.CorrectionSet.from_file(path_to_json)[name]

    if is_correction:
        sf = evaluator.evaluate(events.Pileup.nTrueInt, "nominal")
        sfup, sfdown = None, None

    else:
        sf = np.ones(len(weights._weight))
        sf_nom = evaluator.evaluate(events.Pileup.nTrueInt, "nominal")

        sfup = evaluator.evaluate(events.Pileup.nTrueInt, "up") / sf_nom
        sfdown = evaluator.evaluate(events.Pileup.nTrueInt, "down") / sf_nom

    weights.add(name="Pileup", weight=sf, weightUp=sfup, weightDown=sfdown)

    return weights


def LooseMvaSF(photons, weights, year="2017", is_correction=True, **kwargs):
    """
    LooseMvaSF: correction to the event weight on a per photon level, impacting one of the high importance input variable of the DiphotonBDT, binned in eta and r9.
    for original implementation look at: https://github.com/cms-analysis/flashgg/blob/2677dfea2f0f40980993cade55144636656a8a4f/Systematics/python/flashggDiPhotonSystematics2017_Legacy_cfi.py
    And for presentation on the study: https://indico.cern.ch/event/963617/contributions/4103623/attachments/2141570/3608645/Zee_Validation_UL2017_Update_09112020_Prasant.pdf

    2022: up to this point, it applies the 2017 SF with the new formula for combining the SF for the diphoton candidate.
    """

    # era/year defined as parameter of the function, only 2017 is implemented up to now
    avail_years = ["2016", "2016preVFP", "2016postVFP", "2017", "2018", "2022preEE", "2022postEE"]
    if year not in avail_years:
        print(f"\n WARNING: only LooseMvaSF corrections for the year strings {avail_years} are already implemented! \n Exiting. \n")
        exit()
    elif "2016" in year:
        year = "2016"

    # make this read the 2022 files when available!
    # 2017 file should be renamed with the year in its name...
    json_file = os.path.join(os.path.dirname(__file__), f"JSONs/LooseMvaSF/{year}/LooseMvaSF_{year}.json")
    evaluator = correctionlib.CorrectionSet.from_file(json_file)["LooseMvaSF"]
    if year in ["2016", "2017", "2018"]:
        if is_correction:
            # only calculate correction to nominal weight
            sf_lead = evaluator.evaluate(
                "nominal", photons["pho_lead"].ScEta, photons["pho_lead"].r9
            )
            sf_sublead = evaluator.evaluate(
                "nominal", photons["pho_sublead"].ScEta, photons["pho_sublead"].r9
            )
            sf = sf_lead * sf_sublead

            sfup, sfdown = None, None

        else:
            # only calculate systs
            sf = np.ones(len(weights._weight))
            sf_lead = evaluator.evaluate(
                "nominal", photons["pho_lead"].ScEta, photons["pho_lead"].r9
            )
            sf_sublead = evaluator.evaluate(
                "nominal", photons["pho_sublead"].ScEta, photons["pho_sublead"].r9
            )
            _sf = sf_lead * sf_sublead

            sfup_lead = evaluator.evaluate(
                "up", photons["pho_lead"].ScEta, photons["pho_lead"].r9
            )
            sfup_sublead = evaluator.evaluate(
                "up", photons["pho_sublead"].ScEta, photons["pho_sublead"].r9
            )
            sfup = sfup_lead * sfup_sublead / _sf

            sfdown_lead = evaluator.evaluate(
                "down", photons["pho_lead"].ScEta, photons["pho_lead"].r9
            )
            sfdown_sublead = evaluator.evaluate(
                "down", photons["pho_sublead"].ScEta, photons["pho_sublead"].r9
            )
            sfdown = sfdown_lead * sfdown_sublead / _sf

    elif "2022" in year:

        if is_correction:
            # only calculate correction to nominal weight
            # ToDo: include pT!!!
            sf_lead_p_lead = evaluator.evaluate(
                "nominal", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9  # photons["pho_lead"].pt
            )
            sf_lead_p_sublead = evaluator.evaluate(
                "nominal", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9  # photons["pho_sublead"].pt
            )
            sf_sublead_p_lead = evaluator.evaluate(
                "nominal", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9  # photons["pho_lead"].pt
            )
            sf_sublead_p_sublead = evaluator.evaluate(
                "nominal", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9  # photons["pho_sublead"].pt
            )
            sf = sf_lead_p_lead * sf_sublead_p_sublead + sf_lead_p_sublead * sf_sublead_p_lead - sf_lead_p_lead * sf_lead_p_sublead

            sfup, sfdown = None, None

        else:
            # only calculate systs
            sf = np.ones(len(weights._weight))

            # get nominal SF to divide it out
            sf_lead_p_lead = evaluator.evaluate(
                "nominal", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9  # photons["pho_lead"].pt
            )
            sf_lead_p_sublead = evaluator.evaluate(
                "nominal", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9  # photons["pho_sublead"].pt
            )
            sf_sublead_p_lead = evaluator.evaluate(
                "nominal", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9  # photons["pho_lead"].pt
            )
            sf_sublead_p_sublead = evaluator.evaluate(
                "nominal", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9  # photons["pho_sublead"].pt
            )
            _sf = sf_lead_p_lead * sf_sublead_p_sublead + sf_lead_p_sublead * sf_sublead_p_lead - sf_lead_p_lead * sf_lead_p_sublead

            # up SF
            sfup_lead_p_lead = evaluator.evaluate(
                "up", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9  # photons["pho_lead"].pt
            )
            sfup_lead_p_sublead = evaluator.evaluate(
                "up", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9  # photons["pho_sublead"].pt
            )
            sfup_sublead_p_lead = evaluator.evaluate(
                "up", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9  # photons["pho_lead"].pt
            )
            sfup_sublead_p_sublead = evaluator.evaluate(
                "up", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9  # photons["pho_sublead"].pt
            )
            sfup = (sfup_lead_p_lead * sfup_sublead_p_sublead + sfup_lead_p_sublead * sfup_sublead_p_lead - sfup_lead_p_lead * sfup_lead_p_sublead) / _sf

            # down SF
            sfdown_lead_p_lead = evaluator.evaluate(
                "down", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9  # photons["pho_lead"].pt
            )
            sfdown_lead_p_sublead = evaluator.evaluate(
                "down", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9  # photons["pho_sublead"].pt
            )
            sfdown_sublead_p_lead = evaluator.evaluate(
                "down", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9  # photons["pho_lead"].pt
            )
            sfdown_sublead_p_sublead = evaluator.evaluate(
                "down", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9  # photons["pho_sublead"].pt
            )
            sfdown = (sfdown_lead_p_lead * sfdown_sublead_p_sublead + sfdown_lead_p_sublead * sfdown_sublead_p_lead - sfdown_lead_p_lead * sfdown_lead_p_sublead) / _sf

    weights.add(name="LooseMvaSF", weight=sf, weightUp=sfup, weightDown=sfdown)

    return weights


def ElectronVetoSF(photons, weights, year="2017", is_correction=True, **kwargs):
    """
    ElectronVetoSF: correction to the event weight on a per photon level, Conversion safe veto efficiency with event counting method: To check if the FSR photons are passing the e-veto or not.
    binned in abs(SCeta) and r9.
    for original implementation look at: https://github.com/cms-analysis/flashgg/blob/2677dfea2f0f40980993cade55144636656a8a4f/Systematics/python/flashggDiPhotonSystematics2017_Legacy_cfi.py
    And for presentation on the study: https://indico.cern.ch/event/961164/contributions/4089584/attachments/2135019/3596299/Zmmg_UL2017%20With%20CorrMC_Hgg%20(02.11.2020).pdf
    """

    # era/year defined as parameter of the function
    avail_years = ["2016", "2016preVFP", "2016postVFP", "2017", "2018", "2022preEE", "2022postEE"]
    if year not in avail_years:
        print(f"\n WARNING: only eVetoSF corrections for the year strings {avail_years} are already implemented! \n Exiting. \n")
        exit()
    elif "2016" in year:
        year = "2016"

    if year in ["2016", "2017", "2018"]:
        # 2017 file should be renamed with the year in its name...
        json_file = os.path.join(os.path.dirname(__file__), f"JSONs/ElectronVetoSF/{year}/eVetoSF_{year}.json")
        evaluator = correctionlib.CorrectionSet.from_file(json_file)["ElectronVetoSF"]
        if is_correction:
            # only calculate correction to nominal weight
            sf_lead = evaluator.evaluate(
                "nominal", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9
            )
            sf_sublead = evaluator.evaluate(
                "nominal", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9
            )
            sf = sf_lead * sf_sublead

            sfup, sfdown = None, None

        else:
            # only calculate systs
            sf = np.ones(len(weights._weight))
            sf_lead = evaluator.evaluate(
                "nominal", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9
            )
            sf_sublead = evaluator.evaluate(
                "nominal", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9
            )
            _sf = sf_lead * sf_sublead

            sfup_lead = evaluator.evaluate(
                "up", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9
            )
            sfup_sublead = evaluator.evaluate(
                "up", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9
            )
            sfup = sfup_lead * sfup_sublead / _sf

            sfdown_lead = evaluator.evaluate(
                "down", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9
            )
            sfdown_sublead = evaluator.evaluate(
                "down", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9
            )
            sfdown = sfdown_lead * sfdown_sublead / _sf

    elif "2022" in year:
        # presentation of SF: https://indico.cern.ch/event/1360961/#173-run-3-electron-veto-sfs
        if year == "2022preEE":
            json_file = os.path.join(os.path.dirname(__file__), "JSONs/ElectronVetoSF/2022/preEE_CSEV_SFcorrections.json")
        if year == "2022postEE":
            json_file = os.path.join(os.path.dirname(__file__), "JSONs/ElectronVetoSF/2022/postEE_CSEV_SFcorrections.json")
        evaluator = correctionlib.CorrectionSet.from_file(json_file)["CSEV_SFs"]

        if is_correction:
            # only calculate correction to nominal weight
            sf_lead = evaluator.evaluate(
                abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9, "nominal"
            )
            sf_sublead = evaluator.evaluate(
                abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9, "nominal"
            )
            sf = sf_lead * sf_sublead

            sfup, sfdown = None, None

        else:
            # only calculate systs
            sf = np.ones(len(weights._weight))
            sf_lead = evaluator.evaluate(
                abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9, "nominal"
            )
            sf_sublead = evaluator.evaluate(
                abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9, "nominal"
            )
            _sf = sf_lead * sf_sublead

            unc_lead = evaluator.evaluate(
                abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9, "uncertainty"
            )
            unc_sublead = evaluator.evaluate(
                abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9, "uncertainty"
            )

            sfup = (sf_lead + unc_lead) * (sf_sublead + unc_sublead) / _sf
            sfdown = (sf_lead - unc_lead) * (sf_sublead - unc_sublead) / _sf

    weights.add(name="ElectronVetoSF", weight=sf, weightUp=sfup, weightDown=sfdown)

    return weights


def PreselSF(photons, weights, year="2017", is_correction=True, **kwargs):
    """
    Preselection SF: correction to the event weight on a per photon level for UL2017. Dt:17/11/2020
    Binned in abs(SCeta) and r9.
    For original implementation look at: https://github.com/cms-analysis/flashgg/blob/2677dfea2f0f40980993cade55144636656a8a4f/Systematics/python/flashggDiPhotonSystematics2017_Legacy_cfi.py
    Link to the Presentation: https://indico.cern.ch/event/963617/contributions/4103623/attachments/2141570/3608645/Zee_Validation_UL2017_Update_09112020_Prasant.pdf
    """

    # era/year defined as parameter of the function
    avail_years = ["2016", "2016preVFP", "2016postVFP", "2017", "2018", "2022preEE", "2022postEE"]
    if year not in avail_years:
        print(f"\n WARNING: only PreselSF corrections for the year strings {avail_years} are already implemented! \n Exiting. \n")
        exit()
    elif "2016" in year:
        year = "2016"

    if year in ["2016", "2017", "2018"]:
        json_file = os.path.join(os.path.dirname(__file__), f"JSONs/Preselection/{year}/PreselSF_{year}.json")
    elif year == "2022preEE":
        json_file = os.path.join(os.path.dirname(__file__), "JSONs/Preselection/2022/Preselection_2022PreEE.json")
    elif year == "2022postEE":
        json_file = os.path.join(os.path.dirname(__file__), "JSONs/Preselection/2022/Preselection_2022PostEE.json")

    if year in ["2016", "2017", "2018"]:
        evaluator = correctionlib.CorrectionSet.from_file(json_file)["PreselSF"]
    elif "2022" in year:
        evaluator = correctionlib.CorrectionSet.from_file(json_file)["Preselection_SF"]

    if year in ["2016", "2017", "2018"]:
        if is_correction:
            # only calculate correction to nominal weight
            sf_lead = evaluator.evaluate(
                "nominal", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9
            )
            sf_sublead = evaluator.evaluate(
                "nominal", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9
            )
            sf = sf_lead * sf_sublead

            sfup, sfdown = None, None

        else:
            # only calculate systs
            sf = np.ones(len(weights._weight))
            sf_lead = evaluator.evaluate(
                "nominal", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9
            )
            sf_sublead = evaluator.evaluate(
                "nominal", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9
            )
            _sf = sf_lead * sf_sublead

            sfup_lead = evaluator.evaluate(
                "up", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9
            )
            sfup_sublead = evaluator.evaluate(
                "up", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9
            )
            sfup = sfup_lead * sfup_sublead / _sf

            sfdown_lead = evaluator.evaluate(
                "down", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9
            )
            sfdown_sublead = evaluator.evaluate(
                "down", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9
            )
            sfdown = sfdown_lead * sfdown_sublead / _sf

    # In principle, we should use the fully correct formula https://indico.cern.ch/event/1360948/contributions/5783762/attachments/2788516/4870824/24_02_02_HIG-23-014_PreAppPres.pdf#page=7
    # However, if the SF is pt-binned, the approximation of the multiplication of the two SFs is fully exact
    elif "2022" in year:
        if is_correction:
            # only calculate correction to nominal weight
            sf_lead = evaluator.evaluate(
                abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9, photons["pho_lead"].pt, "nominal"
            )
            sf_sublead = evaluator.evaluate(
                abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9, photons["pho_sublead"].pt, "nominal"
            )
            sf = sf_lead * sf_sublead

            sfup, sfdown = None, None

        else:
            # only calculate systs

            # Slightly different calculation compared to 2017
            # In the 2022 JSONs, the delta is saved as the uncertainty, not the up/down variations of (SF+-delta) themselves
            # Note that the uncertainty is assumed to be symmetric

            sf = np.ones(len(weights._weight))
            sf_lead = evaluator.evaluate(
                abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9, photons["pho_lead"].pt, "nominal"
            )
            sf_sublead = evaluator.evaluate(
                abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9, photons["pho_sublead"].pt, "nominal"
            )
            _sf = sf_lead * sf_sublead

            sf_unc_lead = evaluator.evaluate(
                abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9, photons["pho_lead"].pt, "uncertainty"
            )
            sf_unc_sublead = evaluator.evaluate(
                abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9, photons["pho_sublead"].pt, "uncertainty"
            )

            sfup = (sf_lead + sf_unc_lead) * (sf_sublead + sf_unc_sublead) / _sf

            sfdown = (sf_lead - sf_unc_lead) * (sf_sublead - sf_unc_sublead) / _sf

    weights.add(name="PreselSF", weight=sf, weightUp=sfup, weightDown=sfdown)

    return weights


def TriggerSF(photons, weights, year="2017", is_correction=True, **kwargs):
    """
    Trigger SF: for full 2017 legacy  B-F dataset. Trigger scale factors for use without HLT applied in MC
    Binned in abs(SCeta), r9 and pt.
    For original implementation look at: https://github.com/cms-analysis/flashgg/blob/2677dfea2f0f40980993cade55144636656a8a4f/Systematics/python/flashggDiPhotonSystematics2017_Legacy_cfi.py
    """

    # era/year defined as parameter of the function
    avail_years = ["2016", "2016preVFP", "2016postVFP", "2017", "2018", "2022preEE", "2022postEE"]
    if year not in avail_years:
        print(f"\n WARNING: only TriggerSF corrections for the year strings {avail_years} are already implemented! \n Exiting. \n")
        exit()
    elif "2016" in year:
        year = "2016"

    jsonpog_file_lead = os.path.join(os.path.dirname(__file__), f"JSONs/TriggerSF/{year}/TriggerSF_lead_{year}.json")
    jsonpog_file_sublead = os.path.join(os.path.dirname(__file__), f"JSONs/TriggerSF/{year}/TriggerSF_sublead_{year}.json")
    evaluator_lead = correctionlib.CorrectionSet.from_file(jsonpog_file_lead)["TriggerSF"]
    evaluator_sublead = correctionlib.CorrectionSet.from_file(jsonpog_file_sublead)["TriggerSF"]

    if year in ["2016", "2017", "2018"]:
        if is_correction:
            # only calculate correction to nominal weight
            sf_lead = evaluator_lead.evaluate(
                "nominal", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9, photons["pho_lead"].pt
            )
            sf_sublead = evaluator_sublead.evaluate(
                "nominal", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9, photons["pho_sublead"].pt
            )
            sf = sf_lead * sf_sublead

            sfup, sfdown = None, None

        else:
            # only calculate systs
            sf = np.ones(len(weights._weight))
            sf_lead = evaluator_lead.evaluate(
                "nominal", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9, photons["pho_lead"].pt
            )
            sf_sublead = evaluator_sublead.evaluate(
                "nominal", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9, photons["pho_sublead"].pt
            )
            _sf = sf_lead * sf_sublead

            sfup_lead = evaluator_lead.evaluate(
                "up", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9, photons["pho_lead"].pt
            )
            sfup_sublead = evaluator_sublead.evaluate(
                "up", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9, photons["pho_sublead"].pt
            )
            sfup = sfup_lead * sfup_sublead / _sf

            sfdown_lead = evaluator_lead.evaluate(
                "down", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9, photons["pho_lead"].pt
            )
            sfdown_sublead = evaluator_sublead.evaluate(
                "down", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9, photons["pho_sublead"].pt
            )
            sfdown = sfdown_lead * sfdown_sublead / _sf

    elif "2022" in year:

        if is_correction:
            # only calculate correction to nominal weight
            # ToDo: include pT!!!
            sf_lead_p_lead = evaluator_lead.evaluate(
                "nominal", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9, photons["pho_lead"].pt
            )
            sf_lead_p_sublead = evaluator_lead.evaluate(
                "nominal", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9, photons["pho_sublead"].pt
            )
            sf_sublead_p_lead = evaluator_sublead.evaluate(
                "nominal", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9, photons["pho_lead"].pt
            )
            sf_sublead_p_sublead = evaluator_sublead.evaluate(
                "nominal", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9, photons["pho_sublead"].pt
            )
            sf = sf_lead_p_lead * sf_sublead_p_sublead + sf_lead_p_sublead * sf_sublead_p_lead - sf_lead_p_lead * sf_lead_p_sublead

            sfup, sfdown = None, None

        else:
            # only calculate systs
            sf = np.ones(len(weights._weight))

            # get nominal SF to divide it out
            sf_lead_p_lead = evaluator_lead.evaluate(
                "nominal", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9, photons["pho_lead"].pt
            )
            sf_lead_p_sublead = evaluator_lead.evaluate(
                "nominal", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9, photons["pho_sublead"].pt
            )
            sf_sublead_p_lead = evaluator_sublead.evaluate(
                "nominal", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9, photons["pho_lead"].pt
            )
            sf_sublead_p_sublead = evaluator_sublead.evaluate(
                "nominal", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9, photons["pho_sublead"].pt
            )
            _sf = sf_lead_p_lead * sf_sublead_p_sublead + sf_lead_p_sublead * sf_sublead_p_lead - sf_lead_p_lead * sf_lead_p_sublead

            # up SF
            sfup_lead_p_lead = evaluator_lead.evaluate(
                "up", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9, photons["pho_lead"].pt
            )
            sfup_lead_p_sublead = evaluator_lead.evaluate(
                "up", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9, photons["pho_sublead"].pt
            )
            sfup_sublead_p_lead = evaluator_sublead.evaluate(
                "up", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9, photons["pho_lead"].pt
            )
            sfup_sublead_p_sublead = evaluator_sublead.evaluate(
                "up", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9, photons["pho_sublead"].pt
            )
            sfup = (sfup_lead_p_lead * sfup_sublead_p_sublead + sfup_lead_p_sublead * sfup_sublead_p_lead - sfup_lead_p_lead * sfup_lead_p_sublead) / _sf

            # down SF
            sfdown_lead_p_lead = evaluator_lead.evaluate(
                "down", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9, photons["pho_lead"].pt
            )
            sfdown_lead_p_sublead = evaluator_lead.evaluate(
                "down", abs(photons["pho_lead"].ScEta), photons["pho_lead"].r9, photons["pho_sublead"].pt
            )
            sfdown_sublead_p_lead = evaluator_sublead.evaluate(
                "down", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9, photons["pho_lead"].pt
            )
            sfdown_sublead_p_sublead = evaluator_sublead.evaluate(
                "down", abs(photons["pho_sublead"].ScEta), photons["pho_sublead"].r9, photons["pho_sublead"].pt
            )
            sfdown = (sfdown_lead_p_lead * sfdown_sublead_p_sublead + sfdown_lead_p_sublead * sfdown_sublead_p_lead - sfdown_lead_p_lead * sfdown_lead_p_sublead) / _sf

    weights.add(name="TriggerSF", weight=sf, weightUp=sfup, weightDown=sfdown)

    return weights


def NNLOPS(
    events, dataset_name, weights, is_correction=True, generator="mcatnlo", **kwargs
):
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
            spline_0jet = interp1d(
                nnlops_reweight["0jet"]["pt"], nnlops_reweight["0jet"]["weight"]
            )
            spline_1jet = interp1d(
                nnlops_reweight["1jet"]["pt"], nnlops_reweight["1jet"]["weight"]
            )
            spline_2jet = interp1d(
                nnlops_reweight["2jet"]["pt"], nnlops_reweight["2jet"]["weight"]
            )
            spline_ge3jet = interp1d(
                nnlops_reweight["3jet"]["pt"], nnlops_reweight["3jet"]["weight"]
            )

            # Load truth Higgs pt and njets (pt>30) from events
            higgs_pt = events.HTXS.Higgs_pt
            njets30 = events.HTXS.njets30

            # Extract scale factors from splines and mask for different jet bins
            # Define maximum pt values as interpolated splines only go up so far
            sf = (
                (njets30 == 0) * spline_0jet(np.minimum(np.array(higgs_pt), 125.0))
                + (njets30 == 1) * spline_1jet(np.minimum(np.array(higgs_pt), 625.0))
                + (njets30 == 2) * spline_2jet(np.minimum(np.array(higgs_pt), 800.0))
                + (njets30 >= 3) * spline_ge3jet(np.minimum(np.array(higgs_pt), 925.0))
            )

        else:
            sf = np.ones(len(weights._weight))

    else:
        raise RuntimeError(
            "NNLOPS reweighting is only a flat correction, not a systematic"
        )

    weights.add("NNLOPS", sf, None, None)

    return weights


def AlphaS(photons, events, weights, dataset_name, **kwargs):
    """
    AlphaS weights variations are the last two of the PDF replicas, e.g.,
    https://github.com/cms-sw/cmssw/blob/d37d2797dffc978a78da2fafec3ba480071a0e67/PhysicsTools/NanoAOD/python/genWeightsTable_cfi.py#L10
    https://lhapdfsets.web.cern.ch/current/NNPDF31_nnlo_as_0118_mc_hessian_pdfas/NNPDF31_nnlo_as_0118_mc_hessian_pdfas.info
    """
    systematic = "AlphaS Weight"
    try:
        weights.add(
            name="AlphaS",
            weight=np.ones(len(events)),
            weightUp=events.LHEPdfWeight[:, -1],
            weightDown=events.LHEPdfWeight[:, -2],
        )
    except:
        logger.debug(
            f"No LHEPdf Weights in dataset {dataset_name}, skip systematic: {systematic}"
        )
        return weights

    return weights


def PartonShower(photons, events, weights, dataset_name, **kwargs):
    """
    Parton Shower weights:
    https://github.com/cms-sw/cmssw/blob/caeae4110ddbada1cfdac195404b3c618584e8fb/PhysicsTools/NanoAOD/plugins/GenWeightsTableProducer.cc#L533-L534
    """
    systematic = "PartonShower weight"
    try:
        weights.add(
            name="PS_ISR",
            weight=np.ones(len(events)),
            weightUp=events.PSWeight[:, 0],
            weightDown=events.PSWeight[:, 2],
        )

        weights.add(
            name="PS_FSR",
            weight=np.ones(len(events)),
            weightUp=events.PSWeight[:, 1],
            weightDown=events.PSWeight[:, 3],
        )
    except:
        logger.debug(
            f"No PS Weights in dataset {dataset_name}, skip systematic: {systematic}"
        )
        return weights

    return weights


def cTagSF(events, weights, is_correction=True, year="2017", **kwargs):
    """
    Add c-tagging reshaping SFs as from /https://github.com/higgs-charm/flashgg/blob/dev/cH_UL_Run2_withBDT/Systematics/scripts/applyCTagCorrections.py
    BTV scale factor Wiki: https://btv-wiki.docs.cern.ch/ScaleFactors/
    events must contain jet objects, moreover evaluation of SFs works by calculating the scale factors for all the jets in the event,
    to do this in columnar style the only thing I could think of was to pad the jet collection to the max(n_jets) keep track of the "fake jets" introduced
    by this procedure and fill these position wit 1s before actually setting the weights in the collection. If someone has better ideas I'm open for suggestions
    """

    # era/year defined as parameter of the function, only Run2 is implemented up to now
    avail_years = ["2016preVFP", "2016postVFP", "2017", "2018"]
    if year not in avail_years:
        print(f"\n WARNING: only scale corrections for the year strings {avail_years} are already implemented! \n Exiting. \n")
        exit()

    ctag_systematics = [
        "Extrap",
        "Interp",
        "LHEScaleWeight_muF",
        "LHEScaleWeight_muR",
        "PSWeightFSR",
        "PSWeightISR",
        "PUWeight",
        "Stat",
        "XSec_BRUnc_DYJets_b",
        "XSec_BRUnc_DYJets_c",
        "XSec_BRUnc_WJets_c",
        "jer",
        "jesTotal",
    ]

    ctag_correction_configs = {
        "2016preVFP": {
            "file": os.path.join(
                os.path.dirname(__file__), "JSONs/cTagSF/2016/ctagging_2016preVFP.json.gz"
            ),
            "method": "deepJet_shape",
            "systs": ctag_systematics,
        },
        "2016postVFP": {
            "file": os.path.join(
                os.path.dirname(__file__), "JSONs/cTagSF/2016/ctagging_2016postVFP.json.gz"
            ),
            "method": "deepJet_shape",
            "systs": ctag_systematics,
        },
        "2017": {
            "file": os.path.join(
                os.path.dirname(__file__), "JSONs/cTagSF/2017/ctagging_2017.json.gz"
            ),
            "method": "deepJet_shape",
            "systs": ctag_systematics,
        },
        "2018": {
            "file": os.path.join(
                os.path.dirname(__file__), "JSONs/cTagSF/2018/ctagging_2018.json.gz"
            ),
            "method": "deepJet_shape",
            "systs": ctag_systematics,
        },
    }

    jsonpog_file = os.path.join(
        os.path.dirname(__file__), ctag_correction_configs[year]["file"]
    )
    evaluator = correctionlib.CorrectionSet.from_file(jsonpog_file)[
        ctag_correction_configs[year]["method"]
    ]

    events["n_jets"] = ak.num(events["sel_jets"])
    max_n_jet = max(events["n_jets"])

    dummy_sf = ak.ones_like(events["event"])

    if is_correction:
        # only calculate correction to nominal weight
        # we will append the scale factors relative to all jets to be multiplied
        _sf = []
        # we need a seres of masks to remember where there were no jets
        masks = []
        # to calculate the SFs we have to distinguish for different number of jets
        for i in range(max_n_jet):
            masks.append(events["n_jets"] > i)

            # I select the nth jet column
            nth_jet_hFlav = choose_jet(events["sel_jets"].hFlav, i, 0)
            nth_jet_DeepFlavour_CvsL = choose_jet(
                events["sel_jets"].btagDeepFlav_CvL, i, 0
            )
            nth_jet_DeepFlavour_CvsB = choose_jet(
                events["sel_jets"].btagDeepFlav_CvB, i, 0
            )
            _sf.append(
                evaluator.evaluate(
                    "central",
                    nth_jet_hFlav,
                    nth_jet_DeepFlavour_CvsL,
                    nth_jet_DeepFlavour_CvsB,
                )
            )

            # and fill the places where we had dummies with ones
            _sf[i] = ak.where(
                masks[i],
                _sf[i],
                dummy_sf,
            )

        sfup, sfdown = None, None
        # here we multiply all the sf for different jets in the event
        sf = dummy_sf
        for nth in _sf:
            sf = sf * nth

        sfs_up = []
        sfs_down = []
        for syst in ctag_systematics:
            sfs_up.append(ak.values_astype(dummy_sf, np.float))
            sfs_down.append(ak.values_astype(dummy_sf, np.float))

        weights.add_multivariation(
            name="cTagSF",
            weight=sf,
            modifierNames=ctag_systematics,
            weightsUp=sfs_up,
            weightsDown=sfs_down,
        )

    else:
        # only calculate correction to nominal weight
        # we will append the scale factors relative to all jets to be multiplied
        _sf = []
        # we need a seres of masks to remember where there were no jets
        masks = []
        # to calculate the SFs we have to distinguish for different number of jets
        for i in range(max_n_jet):
            masks.append(events["n_jets"] > i)

            # I select the nth jet column
            nth_jet_hFlav = choose_jet(events["sel_jets"].hFlav, i, 0)
            nth_jet_DeepFlavour_CvsL = choose_jet(
                events["sel_jets"].btagDeepFlav_CvL, i, 0
            )
            nth_jet_DeepFlavour_CvsB = choose_jet(
                events["sel_jets"].btagDeepFlav_CvB, i, 0
            )
            _sf.append(
                evaluator.evaluate(
                    "central",
                    nth_jet_hFlav,
                    nth_jet_DeepFlavour_CvsL,
                    nth_jet_DeepFlavour_CvsB,
                )
            )

            # and fill the places where we had dummies with ones
            _sf[i] = ak.where(
                masks[i],
                _sf[i],
                dummy_sf,
            )

        # here we multiply all the sf for different jets in the event
        sf = dummy_sf
        for nth in _sf:
            sf = sf * nth

        variations = {}
        for syst_name in ctag_correction_configs[year]["systs"]:
            # we will append the scale factors relative to all jets to be multiplied
            _sfup = []
            _sfdown = []
            variations[syst_name] = {}
            for i in range(max_n_jet):
                # I select the nth jet column
                nth_jet_hFlav = choose_jet(events["sel_jets"].hFlav, i, 0)
                nth_jet_DeepFlavour_CvsL = choose_jet(
                    events["sel_jets"].btagDeepFlav_CvL, i, 0
                )
                nth_jet_DeepFlavour_CvsB = choose_jet(
                    events["sel_jets"].btagDeepFlav_CvB, i, 0
                )

                _sfup.append(
                    evaluator.evaluate(
                        "up_" + syst_name,
                        nth_jet_hFlav,
                        nth_jet_DeepFlavour_CvsL,
                        nth_jet_DeepFlavour_CvsB,
                    )
                )

                _sfdown.append(
                    evaluator.evaluate(
                        "down_" + syst_name,
                        nth_jet_hFlav,
                        nth_jet_DeepFlavour_CvsL,
                        nth_jet_DeepFlavour_CvsB,
                    )
                )

                # and fill the places where we had dummies with ones
                _sfup[i] = ak.where(
                    masks[i],
                    _sfup[i],
                    dummy_sf,
                )
                _sfdown[i] = ak.where(
                    masks[i],
                    _sfdown[i],
                    dummy_sf,
                )
            # here we multiply all the sf for different jets in the event
            sfup = dummy_sf
            sfdown = dummy_sf
            for i in range(len(_sf)):
                sfup = sfup * _sfup[i]
                sfdown = sfdown * _sfdown[i]

            variations[syst_name]["up"] = sfup
            variations[syst_name]["down"] = sfdown

        # coffea weights.add_multivariation() wants a list of arrays for the multiple up and down variations
        sfs_up = [variations[syst_name]["up"] / sf for syst_name in ctag_systematics]
        sfs_down = [
            variations[syst_name]["down"] / sf for syst_name in ctag_systematics
        ]

        weights.add_multivariation(
            name="cTagSF",
            weight=dummy_sf,
            modifierNames=ctag_systematics,
            weightsUp=sfs_up,
            weightsDown=sfs_down,
            shift=False,
        )

    return weights


def Zpt(
    events,
    weights,
    logger,
    dataset_name,
    is_correction=True,
    year="2022postEE",
    **kwargs,
):
    """
    Z pt reweighting
    """
    systematic = "Z pt reweighting"

    json_dict = {
        "2016postVFP_UL": os.path.join(
            os.path.dirname(__file__),
            "./JSONs/my_Zpt_reweighting.json.gz",
        ),
        "2016preVFP_UL": os.path.join(
            os.path.dirname(__file__),
            "./JSONs/my_Zpt_reweighting.json.gz",
        ),
        "2017_UL": os.path.join(
            os.path.dirname(__file__),
            "./JSONs/my_Zpt_reweighting.json.gz",
        ),
        "2018_UL": os.path.join(
            os.path.dirname(__file__),
            "./JSONs/my_Zpt_reweighting.json.gz",
        ),
        "2022postEE": os.path.join(
            os.path.dirname(__file__),
            "./JSONs/my_Zpt_reweighting.json.gz",
        ),
        "2023": os.path.join(
            os.path.dirname(__file__),
            "./JSONs/my_Zpt_reweighting.json.gz",
        ),
    }
    key_map = {
        "2016postVFP_UL": "Zpt_reweight",
        "2016preVFP_UL": "Zpt_reweight",
        "2017_UL": "Zpt_reweight",
        "2018_UL": "Zpt_reweight",
        "2022postEE": "Zpt_reweight",
        "2023": "Zpt_reweight",
    }

    # inputs
    input_value = {
        "Zpt": events.mmy_pt,
    }
    cset = correctionlib.CorrectionSet.from_file(json_dict[year])
    # list(cset) # get keys in cset
    sf = cset[key_map[year]]

    logger.debug(f"{systematic}:{key_map[year]}, year: {year} ===> {dataset_name}")
    if is_correction:
        nom = sf.evaluate(input_value["Zpt"])
        weights.add(name="ZptWeight", weight=nom)
    else:
        nom = sf.evaluate(input_value["Zpt"])
        up = sf.evaluate(input_value["Zpt"])
        down = sf.evaluate(input_value["Zpt"])
        weights.add(
            name="ZptWeight",
            weight=ak.ones_like(nom),
            weightUp=up / nom,
            weightDown=down / nom,
        )

    return weights
