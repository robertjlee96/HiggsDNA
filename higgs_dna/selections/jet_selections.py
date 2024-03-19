from higgs_dna.selections.object_selections import delta_r_mask
import awkward
import correctionlib
import os
from coffea.analysis_tools import PackedSelection
from copy import deepcopy
import numpy as np
from correctionlib.highlevel import model_auto, open_auto
import json


def select_jets(
    self,
    jets: awkward.highlevel.Array,
    diphotons: awkward.highlevel.Array,
    muons: awkward.highlevel.Array,
    electrons: awkward.highlevel.Array,
) -> awkward.highlevel.Array:
    pt_cut = jets.pt > self.jet_pt_threshold
    eta_cut = abs(jets.eta) < self.jet_max_eta
    dr_dipho_cut = awkward.ones_like(pt_cut) > 0
    if self.clean_jet_dipho & (awkward.count(diphotons) > 0):
        dr_dipho_cut = delta_r_mask(jets, diphotons, self.jet_dipho_min_dr)

    if (self.clean_jet_pho) & (awkward.count(diphotons) > 0):
        lead = awkward.zip(
            {
                "pt": diphotons.pho_lead.pt,
                "eta": diphotons.pho_lead.eta,
                "phi": diphotons.pho_lead.phi,
                "mass": diphotons.pho_lead.mass,
                "charge": diphotons.pho_lead.charge,
            }
        )
        lead = awkward.with_name(lead, "PtEtaPhiMCandidate")
        sublead = awkward.zip(
            {
                "pt": diphotons.pho_sublead.pt,
                "eta": diphotons.pho_sublead.eta,
                "phi": diphotons.pho_sublead.phi,
                "mass": diphotons.pho_sublead.mass,
                "charge": diphotons.pho_sublead.charge,
            }
        )
        sublead = awkward.with_name(sublead, "PtEtaPhiMCandidate")
        dr_pho_lead_cut = delta_r_mask(jets, lead, self.jet_pho_min_dr)
        dr_pho_sublead_cut = delta_r_mask(jets, sublead, self.jet_pho_min_dr)
    else:
        dr_pho_lead_cut = jets.pt > -1
        dr_pho_sublead_cut = jets.pt > -1

    if (self.clean_jet_ele) & (awkward.count(electrons) > 0):
        dr_electrons_cut = delta_r_mask(jets, electrons, self.jet_ele_min_dr)
    else:
        dr_electrons_cut = jets.pt > -1

    if (self.clean_jet_muo) & (awkward.count(muons) > 0):
        dr_muons_cut = delta_r_mask(jets, muons, self.jet_muo_min_dr)
    else:
        dr_muons_cut = jets.pt > -1

    return (
        (pt_cut)
        & (eta_cut)
        & (dr_dipho_cut)
        & (dr_pho_lead_cut)
        & (dr_pho_sublead_cut)
        & (dr_electrons_cut)
        & (dr_muons_cut)
    )


def select_fatjets(
    self,
    fatjets: awkward.highlevel.Array,
    diphotons: awkward.highlevel.Array,
    muons: awkward.highlevel.Array,
    electrons: awkward.highlevel.Array,
) -> awkward.highlevel.Array:
    # same as select_jets(), but uses fatjet variables
    pt_cut = fatjets.pt > self.fatjet_pt_threshold
    eta_cut = abs(fatjets.eta) < self.fatjet_max_eta
    dr_dipho_cut = awkward.ones_like(pt_cut) > 0
    if self.clean_fatjet_dipho & (awkward.count(diphotons) > 0):
        dr_dipho_cut = delta_r_mask(fatjets, diphotons, self.fatjet_dipho_min_dr)

    if (self.clean_fatjet_pho) & (awkward.count(diphotons) > 0):
        lead = awkward.zip(
            {
                "pt": diphotons.pho_lead.pt,
                "eta": diphotons.pho_lead.eta,
                "phi": diphotons.pho_lead.phi,
                "mass": diphotons.pho_lead.mass,
                "charge": diphotons.pho_lead.charge,
            }
        )
        lead = awkward.with_name(lead, "PtEtaPhiMCandidate")
        sublead = awkward.zip(
            {
                "pt": diphotons.pho_sublead.pt,
                "eta": diphotons.pho_sublead.eta,
                "phi": diphotons.pho_sublead.phi,
                "mass": diphotons.pho_sublead.mass,
                "charge": diphotons.pho_sublead.charge,
            }
        )
        sublead = awkward.with_name(sublead, "PtEtaPhiMCandidate")
        dr_pho_lead_cut = delta_r_mask(fatjets, lead, self.fatjet_pho_min_dr)
        dr_pho_sublead_cut = delta_r_mask(fatjets, sublead, self.fatjet_pho_min_dr)
    else:
        dr_pho_lead_cut = fatjets.pt > -1
        dr_pho_sublead_cut = fatjets.pt > -1

    if (self.clean_fatjet_ele) & (awkward.count(electrons) > 0):
        dr_electrons_cut = delta_r_mask(fatjets, electrons, self.fatjet_ele_min_dr)
    else:
        dr_electrons_cut = fatjets.pt > -1

    if (self.clean_fatjet_muo) & (awkward.count(muons) > 0):
        dr_muons_cut = delta_r_mask(fatjets, muons, self.fatjet_muo_min_dr)
    else:
        dr_muons_cut = fatjets.pt > -1

    return (
        (pt_cut)
        & (eta_cut)
        & (dr_dipho_cut)
        & (dr_pho_lead_cut)
        & (dr_pho_sublead_cut)
        & (dr_electrons_cut)
        & (dr_muons_cut)
    )


def jetvetomap(events, logger, dataset_name, year="2022preEE"):
    """
    Jet veto map
    """
    systematic = "jetvetomap"
    sel_obj = PackedSelection()

    json_dict = {
        "2016preVFP": os.path.join(
            os.path.dirname(__file__),
            "../systematics/JSONs/POG/JME/2016preVFP_UL/jetvetomaps.json.gz",
        ),
        "2016postVFP": os.path.join(
            os.path.dirname(__file__),
            "../systematics/JSONs/POG/JME/2016postVFP_UL/jetvetomaps.json.gz",
        ),
        "2017": os.path.join(
            os.path.dirname(__file__),
            "../systematics/JSONs/POG/JME/2017_UL/jetvetomaps.json.gz",
        ),
        "2018": os.path.join(
            os.path.dirname(__file__),
            "../systematics/JSONs/POG/JME/2018_UL/jetvetomaps.json.gz",
        ),
        "2022postEE": os.path.join(
            os.path.dirname(__file__),
            "../systematics/JSONs/POG/JME/2022_Prompt/jetvetomaps.json.gz",
        ),
        "2022preEE": os.path.join(
            os.path.dirname(__file__),
            "../systematics/JSONs/POG/JME/2022_Prompt/jetvetomaps.json.gz",
        ),
        "2023": os.path.join(
            os.path.dirname(__file__),
            "../systematics/JSONs/POG/JME/2022_Prompt/jetvetomaps.json.gz",
        ),
    }
    key_map = {
        "2016preVFP": "Summer19UL16_V1",
        "2016postVFP": "Summer19UL16_V1",
        "2017": "Summer19UL17_V1",
        "2018": "Summer19UL18_V1",
        "2022preEE": "Winter22Run3_RunCD_V1",
        "2022postEE": "Winter22Run3_RunE_V1",
        "2023": "Winter22Run3_RunCD_V1",
    }

    logger.debug(
        f"[{systematic}] {key_map[year]}, year: {year} to dataset: {dataset_name}"
    )

    # Edge check of input variables. The eta and phi variables don't enable flow
    # https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/JME_2022_Prompt_jetvetomaps.html
    _cset = model_auto(open_auto(json_dict[year]))
    _cset_json = json.loads(_cset.json())
    low_eta, high_eta = (
        _cset_json["corrections"][0]["data"]["content"][0]["value"]["edges"][0][0],
        _cset_json["corrections"][0]["data"]["content"][0]["value"]["edges"][0][-1],
    )
    # phi value must be within [-np.pi,np.pi]. Though values beyond are observed.
    # Might due to the accuracy of nanoaod format. So clip the values to be within the first and last bin centers
    low_phi, high_phi = (
        (
            _cset_json["corrections"][0]["data"]["content"][0]["value"]["edges"][1][0]
            + _cset_json["corrections"][0]["data"]["content"][0]["value"]["edges"][1][1]
        )
        / 2,
        (
            _cset_json["corrections"][0]["data"]["content"][0]["value"]["edges"][1][-1]
            + _cset_json["corrections"][0]["data"]["content"][0]["value"]["edges"][1][
                -2
            ]
        )
        / 2,
    )
    jets_jagged = deepcopy(events.Jet)
    # remove jets out of bin edges
    # https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/JME_2022_Prompt_jetvetomaps.html
    jets_jagged = jets_jagged[
        (jets_jagged.eta >= low_eta) & (jets_jagged.eta < high_eta)
    ]
    count = awkward.num(jets_jagged)
    jets = awkward.flatten(jets_jagged)

    input_dict = {
        "type": "jetvetomap",
        "eta": jets.eta,
        "phi": np.clip(jets.phi, low_phi, high_phi),
    }
    if year != "2022postEE":
        cset = correctionlib.CorrectionSet.from_file(json_dict[year])
        inputs = [input_dict[input.name] for input in cset[key_map[year]].inputs]
        vetomap = cset[key_map[year]].evaluate(*(inputs))
        sel_obj.add("vetomap", np.abs(vetomap) > 0)
    else:
        # ref: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis#From_JME
        # normal jetvetomap should be applied to all 2022
        year = "2022preEE"
        cset = correctionlib.CorrectionSet.from_file(json_dict["2022preEE"])
        inputs = [input_dict[input.name] for input in cset[key_map["2022preEE"]].inputs]
        vetomap = cset[key_map["2022preEE"]].evaluate(*(inputs))

        # consider the EELeak region
        year = "2022postEE"
        cset_eep = correctionlib.CorrectionSet.from_file(json_dict[year])
        input_dict["type"] = "jetvetomap_eep"
        inputs = [input_dict[input.name] for input in cset_eep[key_map[year]].inputs]
        vetomap_eep = cset_eep[key_map[year]].evaluate(*(inputs))
        sel_obj.add(
            "vetomap",
            (np.abs(vetomap) > 0) | ((np.abs(vetomap_eep) > 0) & (jets.pt > 30)),
        )

    sel_veto_jet = sel_obj.all(*(sel_obj.names))
    sel_good_jet = ~awkward.Array(sel_veto_jet)
    logger.info(
        f"[{systematic}] total: {len(sel_good_jet)}, pass: {awkward.sum(sel_good_jet)}"
    )
    sel_good_jet_jagged = awkward.unflatten(sel_good_jet, count)
    events.Jet = jets_jagged[sel_good_jet_jagged]

    return events
