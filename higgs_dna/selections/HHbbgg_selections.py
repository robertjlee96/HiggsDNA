import awkward as ak
import numpy as np
import vector

vector.register_awkward()


def get_HHbbgg(
    self, diphotons: ak.highlevel.Array, dijets: ak.highlevel.Array
) -> ak.highlevel.Array:
    # Script adapted from the Zmumug analysis
    # combine dijet & diphoton
    dijets["charge"] = ak.zeros_like(dijets.pt, dtype=np.int)
    HHbbgg_jagged = ak.cartesian({"diphoton": diphotons, "dijet": dijets}, axis=1)
    # flatten HHbbgg, selection only accept flatten arrays
    count = ak.num(HHbbgg_jagged)
    HHbbgg = ak.flatten(HHbbgg_jagged)

    # diphoton and dijet obj
    diphoton_obj = HHbbgg.diphoton.pho_lead + HHbbgg.diphoton.pho_sublead
    dijet_obj = HHbbgg.dijet.first_jet + HHbbgg.dijet.second_jet

    HHbbgg_obj = HHbbgg.diphoton.pho_lead + HHbbgg.diphoton.pho_sublead + HHbbgg.dijet.first_jet + HHbbgg.dijet.second_jet

    # dress other variables
    HHbbgg["obj_diphoton"] = diphoton_obj
    HHbbgg["pho_lead"] = HHbbgg.diphoton.pho_lead
    HHbbgg["pho_sublead"] = HHbbgg.diphoton.pho_sublead
    HHbbgg["obj_dijet"] = dijet_obj
    HHbbgg["first_jet"] = HHbbgg.dijet.first_jet
    HHbbgg["second_jet"] = HHbbgg.dijet.second_jet
    HHbbgg["obj_HHbbgg"] = HHbbgg_obj
    HHbbgg_jagged = ak.unflatten(HHbbgg, count)

    # get best matched HHbbgg for each event
    best_HHbbgg = ak.firsts(HHbbgg_jagged)

    return best_HHbbgg


def getCosThetaStar_CS(HHbbgg, ebeam):
    """
    cos theta star angle in the Collins Soper frame
    Copied directly from here: https://github.com/ResonantHbbHgg/Selection/blob/master/selection.h#L3367-L3385
    """
    p1 = ak.zip(
        {
            "px": 0,
            "py": 0,
            "pz": ebeam,
            "E": ebeam,
        },
        with_name="Momentum4D",
    )

    p2 = ak.zip(
        {
            "px": 0,
            "py": 0,
            "pz": -ebeam,
            "E": ebeam,
        },
        with_name="Momentum4D",
    )

    diphoton = ak.with_name(HHbbgg.obj_diphoton,"Momentum4D")
    # dijet=ak.with_name(HHbbgg.obj_dijet,"Momentum4D")
    HH = ak.with_name(HHbbgg.obj_HHbbgg,"Momentum4D")

    hhforboost = ak.zip({"px": -HH.px ,"py":-HH.py, "pz":-HH.pz, "E": HH.E})
    hhforboost = ak.with_name(hhforboost,"Momentum4D")

    p1 = p1.boost(hhforboost)
    p2 = p2.boost(hhforboost)
    diphotonBoosted = diphoton.boost(hhforboost)

    CSaxis = (p1 - p2)

    return np.cos(CSaxis.deltaangle(diphotonBoosted))


def getCosThetaStar_gg(HHbbgg):

    Hgg = ak.with_name(HHbbgg.obj_diphoton,"Momentum4D")
    hggforboost = ak.zip({"px": -Hgg.px ,"py":-Hgg.py, "pz":-Hgg.pz, "E": Hgg.E})
    hggforboost = ak.with_name(hggforboost,"Momentum4D")

    photon = ak.zip({"px": HHbbgg.pho_lead.px ,"py":HHbbgg.pho_lead.py, "pz":HHbbgg.pho_lead.pz, "mass": 0})
    Hgg_photon = ak.with_name(photon,"Momentum4D")

    Hgg_photon_boosted = Hgg_photon.boost(hggforboost)

    return Hgg_photon_boosted.costheta


def getCosThetaStar_jj(HHbbgg):

    Hjj = ak.with_name(HHbbgg.obj_dijet,"Momentum4D")
    hjjforboost = ak.zip({"px": -Hjj.px ,"py":-Hjj.py, "pz":-Hjj.pz, "E": Hjj.E})
    hjjforboost = ak.with_name(hjjforboost,"Momentum4D")

    Hjj_jet = ak.with_name(HHbbgg.first_jet,"Momentum4D")

    Hjj_jet_boosted = Hjj_jet.boost(hjjforboost)

    return Hjj_jet_boosted.costheta


def DeltaR(photon, jet):
    pho_obj = ak.with_name(photon,"Momentum4D")
    jet_obj = ak.with_name(jet,"Momentum4D")
    return vector.Spatial.deltaR(pho_obj,jet_obj)
