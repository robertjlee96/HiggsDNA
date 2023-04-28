import numpy as np
import correctionlib
import os


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

def AlphaS(photons, events, weights, logger, dataset, systematic):
    """
    AlphaS weights variations are the last two of the PDF replicas, e.g.,
    https://github.com/cms-sw/cmssw/blob/d37d2797dffc978a78da2fafec3ba480071a0e67/PhysicsTools/NanoAOD/python/genWeightsTable_cfi.py#L10
    https://lhapdfsets.web.cern.ch/current/NNPDF31_nnlo_as_0118_mc_hessian_pdfas/NNPDF31_nnlo_as_0118_mc_hessian_pdfas.info
    """
    try:
        weights.add(
            name="AlphaS",
            weight=np.ones(len(events)),
            weightUp=events.LHEPdfWeight[:,-1],
            weightDown=events.LHEPdfWeight[:,-2]
        )
    except:
        logger.info(
            f"No LHEPdf Weights in dataset {dataset}, skip systematic {systematic}"
        )
        return weights

    return weights

def PartonShower(photons, events, weights, logger, dataset, systematic):
    """
    AlphaS weights variations are the last two of the PDF replicas, e.g.,
    https://github.com/cms-sw/cmssw/blob/d37d2797dffc978a78da2fafec3ba480071a0e67/PhysicsTools/NanoAOD/python/genWeightsTable_cfi.py#L10
    https://lhapdfsets.web.cern.ch/current/NNPDF31_nnlo_as_0118_mc_hessian_pdfas/NNPDF31_nnlo_as_0118_mc_hessian_pdfas.info
    """
    try:
        weights.add(
            name="PS_FSR",
            weight=np.ones(len(events)),
            weightUp=events.PSWeight[:,0],
            weightDown=events.PSWeight[:,2]
        )

        weights.add(
            name="PS_ISR",
            weight=np.ones(len(events)),
            weightUp=events.PSWeight[:,1],
            weightDown=events.PSWeight[:,3]
        )
    except:
        logger.info(
            f"No PS Weights in dataset {dataset}, skip systematic {systematic}"
        )
        return weights

    return weights


