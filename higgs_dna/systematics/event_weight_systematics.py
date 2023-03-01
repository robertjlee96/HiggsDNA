import numpy as np
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
