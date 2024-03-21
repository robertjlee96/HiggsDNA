from higgs_dna.workflows.base import HggBaseProcessor
from higgs_dna.tools.SC_eta import add_photon_SC_eta
from higgs_dna.tools.EELeak_region import veto_EEleak_flag
from higgs_dna.selections.photon_selections import photon_preselection
from higgs_dna.selections.lepton_selections import select_electrons, select_muons
from higgs_dna.selections.jet_selections import select_jets, select_fatjets, jetvetomap
from higgs_dna.selections.lumi_selections import select_lumis
from higgs_dna.utils.dumping_utils import diphoton_ak_array, dump_ak_array, diphoton_list_to_pandas, dump_pandas
from higgs_dna.utils.misc_utils import choose_jet
from higgs_dna.selections.HHbbgg_selections import getCosThetaStar_CS, get_HHbbgg, getCosThetaStar_gg, getCosThetaStar_jj, DeltaR
from higgs_dna.tools.flow_corrections import calculate_flow_corrections

from higgs_dna.tools.mass_decorrelator import decorrelate_mass_resolution
from higgs_dna.tools.truth_info import get_truth_info_dict

# from higgs_dna.utils.dumping_utils import diphoton_list_to_pandas, dump_pandas
from higgs_dna.systematics import object_systematics as available_object_systematics
from higgs_dna.systematics import object_corrections as available_object_corrections
from higgs_dna.systematics import weight_systematics as available_weight_systematics
from higgs_dna.systematics import weight_corrections as available_weight_corrections

import functools
import warnings
from typing import Any, Dict, List, Optional
import awkward
import numpy
import sys
import vector
from coffea.analysis_tools import Weights
from copy import deepcopy

import logging

logger = logging.getLogger(__name__)

vector.register_awkward()


class HHbbggProcessor(HggBaseProcessor):
    def __init__(
        self,
        metaconditions: Dict[str, Any],
        systematics: Dict[str, List[Any]] = None,
        corrections: Dict[str, List[Any]] = None,
        apply_trigger: bool = False,
        output_location: Optional[str] = None,
        taggers: Optional[List[Any]] = None,
        skipCQR: bool = False,
        skipJetVetoMap: bool = False,
        year: Dict[str, List[str]] = None,
        fiducialCuts: str = "classical",
        doDeco: bool = False,
        Smear_sigma_m: bool = False,
        doFlow_corrections: bool = False,
        output_format: str = "parquet"
    ) -> None:
        super().__init__(
            metaconditions,
            systematics=systematics,
            corrections=corrections,
            apply_trigger=apply_trigger,
            output_location=output_location,
            taggers=taggers,
            trigger_group=".*DoubleEG.*",
            analysis="mainAnalysis",
            skipCQR=skipCQR,
            skipJetVetoMap=skipJetVetoMap,
            year=year,
            fiducialCuts=fiducialCuts,
            doDeco=doDeco,
            Smear_sigma_m=Smear_sigma_m,
            doFlow_corrections=doFlow_corrections,
            output_format=output_format
        )
        self.trigger_group = ".*DoubleEG.*"
        self.analysis = "mainAnalysis"

        self.num_fatjets_to_store = 4
        self.num_leptons_to_store = 4

        # fatjet selection cuts, same as jets but for photon dR and pT
        self.fatjet_dipho_min_dr = 0.8
        self.fatjet_pho_min_dr = 0.8
        self.fatjet_ele_min_dr = 0.4
        self.fatjet_muo_min_dr = 0.4
        self.fatjet_pt_threshold = 250
        self.fatjet_max_eta = 4.7

        self.clean_fatjet_dipho = False
        self.clean_fatjet_pho = True
        self.clean_fatjet_ele = False
        self.clean_fatjet_muo = False

        # Redefine some cuts
        self.jet_max_eta = 4.7
        self.num_jets_to_store = 6

    def process_extra(self, events: awkward.Array) -> awkward.Array:
        return events, {}

    def process(self, events: awkward.Array) -> Dict[Any, Any]:
        dataset_name = events.metadata["dataset"]
        filename = events.metadata["filename"]

        # data or monte carlo?
        self.data_kind = "mc" if hasattr(events, "GenPart") else "data"

        # here we start recording possible coffea accumulators
        # most likely histograms, could be counters, arrays, ...
        histos_etc = {}
        histos_etc[dataset_name] = {}
        if self.data_kind == "mc":
            histos_etc[dataset_name]["nTot"] = int(
                awkward.num(events.genWeight, axis=0)
            )
            histos_etc[dataset_name]["nPos"] = int(awkward.sum(events.genWeight > 0))
            histos_etc[dataset_name]["nNeg"] = int(awkward.sum(events.genWeight < 0))
            histos_etc[dataset_name]["nEff"] = int(
                histos_etc[dataset_name]["nPos"] - histos_etc[dataset_name]["nNeg"]
            )
            histos_etc[dataset_name]["genWeightSum"] = float(
                awkward.sum(events.genWeight)
            )
        else:
            histos_etc[dataset_name]["nTot"] = int(len(events))
            histos_etc[dataset_name]["nPos"] = int(histos_etc[dataset_name]["nTot"])
            histos_etc[dataset_name]["nNeg"] = int(0)
            histos_etc[dataset_name]["nEff"] = int(histos_etc[dataset_name]["nTot"])
            histos_etc[dataset_name]["genWeightSum"] = float(len(events))

        # lumi mask
        if self.data_kind == "data":
            try:
                lumimask = select_lumis(self.year[dataset_name][0], events, logger)
                events = events[lumimask]
            except:
                logger.info(
                    f"[ lumimask ] Skip now! Unable to find year info of {dataset_name}"
                )
        # apply jetvetomap
        if not self.skipJetVetoMap:
            events = jetvetomap(
                events, logger, dataset_name, year=self.year[dataset_name][0]
            )
        # metadata array to append to higgsdna output
        metadata = {}

        if self.data_kind == "mc":
            # Add sum of gen weights before selection for normalisation in postprocessing
            metadata["sum_genw_presel"] = str(awkward.sum(events.genWeight))
        else:
            metadata["sum_genw_presel"] = "Data"

        # apply filters and triggers
        events = self.apply_filters_and_triggers(events)

        # we need ScEta for corrections and systematics, which is not present in NanoAODv11 but can be calculated using PV
        events.Photon = add_photon_SC_eta(events.Photon, events.PV)

        # add veto EE leak branch for photons, could also be used for electrons
        if (
            self.year[dataset_name][0] == "2022EE"
            or self.year[dataset_name][0] == "2022postEE"
        ):
            events.Photon = veto_EEleak_flag(self, events.Photon)

        # read which systematics and corrections to process
        try:
            correction_names = self.corrections[dataset_name]
        except KeyError:
            correction_names = []
        try:
            systematic_names = self.systematics[dataset_name]
        except KeyError:
            systematic_names = []

        # If --Smear_sigma_m == True and no Smearing correction in .json for MC throws an error, since the pt scpectrum need to be smeared in order to properly calculate the smeared sigma_m_m
        if self.data_kind == "mc" and self.Smear_sigma_m and 'Smearing' not in correction_names:
            warnings.warn("Smearing should be specified in the corrections field in .json in order to smear the mass!")
            sys.exit(0)

        # Since now we are applying Smearing term to the sigma_m_over_m i added this portion of code
        # specially for the estimation of smearing terms for the data events [data pt/energy] are not smeared!
        if self.data_kind == "data" and self.Smear_sigma_m:
            correction_name = 'Smearing'

            logger.info(
                f"\nApplying correction {correction_name} to dataset {dataset_name}\n"
            )
            varying_function = available_object_corrections[correction_name]
            events = varying_function(events=events, year=self.year[dataset_name][0])

        for correction_name in correction_names:
            if correction_name in available_object_corrections.keys():
                logger.info(
                    f"Applying correction {correction_name} to dataset {dataset_name}"
                )
                varying_function = available_object_corrections[correction_name]
                events = varying_function(events=events, year=self.year[dataset_name][0])
            elif correction_name in available_weight_corrections:
                # event weight corrections will be applied after photon preselection / application of further taggers
                continue
            else:
                # may want to throw an error instead, needs to be discussed
                warnings.warn(f"Could not process correction {correction_name}.")
                continue

        original_photons = events.Photon
        original_jets = events.Jet
        # systematic object variations
        for systematic_name in systematic_names:
            if systematic_name in available_object_systematics.keys():
                systematic_dct = available_object_systematics[systematic_name]
                if systematic_dct["object"] == "Photon":
                    logger.info(
                        f"Adding systematic {systematic_name} to photons collection of dataset {dataset_name}"
                    )
                    original_photons.add_systematic(
                        # passing the arguments here explicitly since I want to pass the events to the varying function. If there is a more elegant / flexible way, just change it!
                        name=systematic_name,
                        kind=systematic_dct["args"]["kind"],
                        what=systematic_dct["args"]["what"],
                        varying_function=functools.partial(
                            systematic_dct["args"]["varying_function"],
                            events=events,
                            year=self.year[dataset_name][0],
                        )
                        # name=systematic_name, **systematic_dct["args"]
                    )
                elif systematic_dct["object"] == "Jet":
                    logger.info(
                        f"Adding systematic {systematic_name} to jets collection of dataset {dataset_name}"
                    )
                    original_jets.add_systematic(
                        # passing the arguments here explicitly since I want to pass the events to the varying function. If there is a more elegant / flexible way, just change it!
                        name=systematic_name,
                        kind=systematic_dct["args"]["kind"],
                        what=systematic_dct["args"]["what"],
                        varying_function=functools.partial(
                            systematic_dct["args"]["varying_function"], events=events
                        )
                        # name=systematic_name, **systematic_dct["args"]
                    )
                # to be implemented for other objects here
            elif systematic_name in available_weight_systematics:
                # event weight systematics will be applied after photon preselection / application of further taggers
                continue
            else:
                # may want to throw an error instead, needs to be discussed
                warnings.warn(
                    f"Could not process systematic variation {systematic_name}."
                )
                continue

        # Computing the normalizinf flow correction
        if self.data_kind == "mc" and self.doFlow_corrections:

            # Applyting the Flow corrections to all photons before pre-selection
            counts = awkward.num(original_photons)
            corrected_inputs,var_list = calculate_flow_corrections(original_photons, events, self.meta["flashggPhotons"]["flow_inputs"], self.meta["flashggPhotons"]["Isolation_transform_order"], year=self.year[dataset_name][0])

            # Store the raw nanoAOD value and update photon ID MVA value for preselection
            original_photons["mvaID_run3"] = awkward.unflatten(self.add_photonid_mva_run3(original_photons, events), counts)
            original_photons["mvaID_nano"] = original_photons["mvaID"]

            # Store the raw values of the inputs and update the input values with the corrections since some variables used in the preselection
            for i in range(len(var_list)):
                original_photons["raw_" + str(var_list[i])] = original_photons[str(var_list[i])]
                original_photons[str(var_list[i])] = awkward.unflatten(corrected_inputs[:,i] , counts)

            original_photons["mvaID"] = awkward.unflatten(self.add_photonid_mva_run3(original_photons, events), counts)

        # Applying systematic variations
        photons_dct = {}
        photons_dct["nominal"] = original_photons
        logger.debug(original_photons.systematics.fields)
        for systematic in original_photons.systematics.fields:
            for variation in original_photons.systematics[systematic].fields:
                # deepcopy to allow for independent calculations on photon variables with CQR
                photons_dct[f"{systematic}_{variation}"] = deepcopy(
                    original_photons.systematics[systematic][variation]
                )

        jets_dct = {}
        jets_dct["nominal"] = original_jets
        logger.debug(original_jets.systematics.fields)
        for systematic in original_jets.systematics.fields:
            for variation in original_jets.systematics[systematic].fields:
                # deepcopy to allow for independent calculations on photon variables with CQR
                jets_dct[f"{systematic}_{variation}"] = original_jets.systematics[
                    systematic
                ][variation]

        for variation, photons in photons_dct.items():
            for jet_variation, Jets in jets_dct.items():
                # make sure no duplicate executions
                if variation == "nominal" or jet_variation == "nominal":
                    if variation != "nominal" and jet_variation != "nominal":
                        continue
                    do_variation = "nominal"
                    if not (variation == "nominal" and jet_variation == "nominal"):
                        do_variation = (
                            variation if variation != "nominal" else jet_variation
                        )
                    logger.debug("Variation: {}".format(do_variation))
                    if self.chained_quantile is not None:
                        photons = self.chained_quantile.apply(photons, events)
                    # recompute photonid_mva on the fly
                    if self.photonid_mva_EB and self.photonid_mva_EE:
                        photons = self.add_photonid_mva(photons, events)

                    # photon preselection
                    photons = photon_preselection(self, photons, events, year=self.year[dataset_name][0])
                    # sort photons in each event descending in pt
                    # make descending-pt combinations of photons
                    photons = photons[awkward.argsort(photons.pt, ascending=False)]
                    photons["charge"] = awkward.zeros_like(
                        photons.pt
                    )  # added this because charge is not a property of photons in nanoAOD v11. We just assume every photon has charge zero...
                    diphotons = awkward.combinations(
                        photons, 2, fields=["pho_lead", "pho_sublead"]
                    )
                    # the remaining cut is to select the leading photons
                    # the previous sort assures the order
                    diphotons = diphotons[
                        diphotons["pho_lead"].pt > self.min_pt_lead_photon
                    ]

                    # now turn the diphotons into candidates with four momenta and such
                    diphoton_4mom = diphotons["pho_lead"] + diphotons["pho_sublead"]
                    diphotons["pt"] = diphoton_4mom.pt
                    diphotons["eta"] = diphoton_4mom.eta
                    diphotons["phi"] = diphoton_4mom.phi
                    diphotons["mass"] = diphoton_4mom.mass
                    diphotons["charge"] = diphoton_4mom.charge
                    diphotons = awkward.with_name(diphotons, "PtEtaPhiMCandidate")

                    # sort diphotons by pT
                    diphotons = diphotons[
                        awkward.argsort(diphotons.pt, ascending=False)
                    ]

                    # Determine if event passes fiducial Hgg cuts at detector-level
                    if self.fiducialCuts == 'classical':
                        fid_det_passed = (diphotons.pho_lead.pt / diphotons.mass > 1 / 3) & (diphotons.pho_sublead.pt / diphotons.mass > 1 / 4) & (diphotons.pho_lead.pfRelIso03_all_quadratic * diphotons.pho_lead.pt < 10) & ((diphotons.pho_sublead.pfRelIso03_all_quadratic * diphotons.pho_sublead.pt) < 10) & (numpy.abs(diphotons.pho_lead.eta) < 2.5) & (numpy.abs(diphotons.pho_sublead.eta) < 2.5)
                    elif self.fiducialCuts == 'geometric':
                        fid_det_passed = (numpy.sqrt(diphotons.pho_lead.pt * diphotons.pho_sublead.pt) / diphotons.mass > 1 / 3) & (diphotons.pho_sublead.pt / diphotons.mass > 1 / 4) & (diphotons.pho_lead.pfRelIso03_all_quadratic * diphotons.pho_lead.pt < 10) & (diphotons.pho_sublead.pfRelIso03_all_quadratic * diphotons.pho_sublead.pt < 10) & (numpy.abs(diphotons.pho_lead.eta) < 2.5) & (numpy.abs(diphotons.pho_sublead.eta) < 2.5)
                    elif self.fiducialCuts == 'none':
                        fid_det_passed = diphotons.pho_lead.pt > -10  # This is a very dummy way but I do not know how to make a true array of outer shape of diphotons
                    else:
                        warnings.warn("You chose %s the fiducialCuts mode, but this is currently not supported. You should check your settings. For this run, no fiducial selection at detector level is applied." % self.fiducialCuts)
                        fid_det_passed = diphotons.pho_lead.pt > -10

                    diphotons = diphotons[fid_det_passed]

                    # baseline modifications to diphotons
                    if self.diphoton_mva is not None:
                        diphotons = self.add_diphoton_mva(diphotons, events)

                    # workflow specific processing
                    events, process_extra = self.process_extra(events)
                    histos_etc.update(process_extra)

                    # jet_variables
                    jets = awkward.zip(
                        {
                            "pt": Jets.pt,
                            "eta": Jets.eta,
                            "phi": Jets.phi,
                            "mass": Jets.mass,
                            "charge": awkward.zeros_like(
                                Jets.pt
                            ),  # added this because jet charge is not a property of photons in nanoAOD v11. We just need the charge to build jet collection.
                            "hFlav": Jets.hadronFlavour
                            if self.data_kind == "mc"
                            else awkward.zeros_like(Jets.pt),
                            "btagDeepFlav_B": Jets.btagDeepFlavB,
                            "btagDeepFlav_CvB": Jets.btagDeepFlavCvB,
                            "btagDeepFlav_CvL": Jets.btagDeepFlavCvL,
                            "btagDeepFlav_QG": Jets.btagDeepFlavQG,
                            "btagPNetB": Jets.btagPNetB,
                            "PNetRegPtRawCorr": Jets.PNetRegPtRawCorr,
                            "PNetRegPtRawCorrNeutrino": Jets.PNetRegPtRawCorrNeutrino,
                            "PNetRegPtRawRes": Jets.PNetRegPtRawRes,
                            "btagRobustParTAK4B": Jets.btagRobustParTAK4B,
                        }
                    )
                    jets = awkward.with_name(jets, "PtEtaPhiMCandidate")

                    electrons = awkward.zip(
                        {
                            "pt": events.Electron.pt,
                            "eta": events.Electron.eta,
                            "phi": events.Electron.phi,
                            "mass": events.Electron.mass,
                            "charge": events.Electron.charge,
                            "mvaIso_WP90": events.Electron.mvaIso_WP90,
                            "mvaIso_WP80": events.Electron.mvaIso_WP80,
                        }
                    )
                    electrons = awkward.with_name(electrons, "PtEtaPhiMCandidate")

                    muons = awkward.zip(
                        {
                            "pt": events.Muon.pt,
                            "eta": events.Muon.eta,
                            "phi": events.Muon.phi,
                            "mass": events.Muon.mass,
                            "charge": events.Muon.charge,
                            "tightId": events.Muon.tightId,
                            "mediumId": events.Muon.mediumId,
                            "looseId": events.Muon.looseId,
                            "isGlobal": events.Muon.isGlobal,
                        }
                    )
                    muons = awkward.with_name(muons, "PtEtaPhiMCandidate")

                    # create PuppiMET objects
                    puppiMET = events.PuppiMET
                    puppiMET = awkward.with_name(puppiMET, "PtEtaPhiMCandidate")

                    # FatJet variables
                    fatjets = events.FatJet
                    fatjets["charge"] = awkward.zeros_like(fatjets.pt)
                    fatjets = awkward.with_name(fatjets, "PtEtaPhiMCandidate")

                    # SubJet variables
                    subjets = events.SubJet
                    subjets["charge"] = awkward.zeros_like(subjets.pt)
                    subjets = awkward.with_name(subjets, "PtEtaPhiMCandidate")

                    # GenJetAK8 variables
                    if self.data_kind == "mc":
                        genjetsAK8 = events.GenJetAK8
                        genjetsAK8["charge"] = awkward.zeros_like(genjetsAK8.pt)
                        genjetsAK8 = awkward.with_name(genjetsAK8, "PtEtaPhiMCandidate")

                    # lepton cleaning
                    sel_electrons = electrons[
                        select_electrons(self, electrons, diphotons)
                    ]

                    sel_muons = muons[select_muons(self, muons, diphotons)]

                    # jet selection and pt ordering
                    jets = jets[
                        select_jets(self, jets, diphotons, sel_muons, sel_electrons)
                    ]
                    jets = jets[awkward.argsort(jets.pt, ascending=False)]

                    # fatjet selection and pt ordering
                    fatjets = fatjets[select_fatjets(self, fatjets, diphotons, sel_muons, sel_electrons)]  # For now, having the same preselection as jet. Can be changed later
                    fatjets = fatjets[awkward.argsort(fatjets.particleNetWithMass_HbbvsQCD, ascending=False)]

                    # adding selected jets to events to be used in ctagging SF calculation
                    events["sel_jets"] = jets
                    n_jets = awkward.num(jets)

                    diphotons["n_jets"] = n_jets

                    ## --------- Beginning of the part added for the HHTobbgg analysis -------
                    # Store  6 jets with their infos -> Taken from the top workflow. This part of the code was orignally written by Florain Mausolf
                    jets_properties = ["pt", "eta", "phi", "mass", "charge", "btagDeepFlav_B", "btagPNetB", "PNetRegPtRawCorr", "PNetRegPtRawCorrNeutrino", "PNetRegPtRawRes", "btagRobustParTAK4B"]
                    for i in range(self.num_jets_to_store):  # Number of jets to select
                        for prop in jets_properties:
                            key = f"jet{i+1}_{prop}"
                            # Retrieve the value using the choose_jet function
                            value = choose_jet(getattr(jets, prop), i, -999.0)
                            # Store the value in the diphotons dictionary
                            diphotons[key] = value

                    # Add the truth information
                    if self.data_kind == "mc":
                        param_values = get_truth_info_dict(filename)
                        for key in param_values.keys():
                            diphotons[key] = param_values[key]

                    # Creatiion a dijet
                    dijets = awkward.combinations(
                        jets, 2, fields=("first_jet", "second_jet")
                    )

                    # HHbbgg :  now turn the dijets into candidates with four momenta and such
                    dijets_4mom = dijets["first_jet"] + dijets["second_jet"]
                    dijets["pt"] = dijets_4mom.pt
                    dijets["eta"] = dijets_4mom.eta
                    dijets["phi"] = dijets_4mom.phi
                    dijets["mass"] = dijets_4mom.mass
                    dijets["charge"] = dijets_4mom.charge
                    dijets["btagPNetB_sum"] = dijets["first_jet"].btagPNetB + dijets["second_jet"].btagPNetB
                    dijets = awkward.with_name(dijets, "PtEtaPhiMCandidate")

                    # Selection on the dijet
                    dijets = dijets[dijets.mass < 190]
                    dijets = dijets[dijets.mass > 70]
                    dijets = dijets[dijets.btagPNetB_sum > 0]
                    dijets = dijets[(numpy.abs(dijets["first_jet"].eta) < 2.5) & (numpy.abs(dijets["second_jet"].eta) < 2.5)]
                    dijets = dijets[awkward.argsort(dijets.btagPNetB_sum, ascending=False)]

                    lead_bjet_pt = choose_jet(dijets["first_jet"].pt, 0, -999.0)
                    lead_bjet_eta = choose_jet(dijets["first_jet"].eta, 0, -999.0)
                    lead_bjet_phi = choose_jet(dijets["first_jet"].phi, 0, -999.0)
                    lead_bjet_mass = choose_jet(dijets["first_jet"].mass, 0, -999.0)
                    lead_bjet_charge = choose_jet(dijets["first_jet"].charge, 0, -999.0)
                    lead_bjet_btagPNetB = choose_jet(dijets["first_jet"].btagPNetB, 0, -999.0)
                    lead_bjet_PNetRegPtRawCorr = choose_jet(dijets["first_jet"].PNetRegPtRawCorr, 0, -999.0)
                    lead_bjet_PNetRegPtRawCorrNeutrino = choose_jet(dijets["first_jet"].PNetRegPtRawCorrNeutrino, 0, -999.0)
                    lead_bjet_PNetRegPtRawRes = choose_jet(dijets["first_jet"].PNetRegPtRawRes, 0, -999.0)

                    sublead_bjet_pt = choose_jet(dijets["second_jet"].pt, 0, -999.0)
                    sublead_bjet_eta = choose_jet(dijets["second_jet"].eta, 0, -999.0)
                    sublead_bjet_phi = choose_jet(dijets["second_jet"].phi, 0, -999.0)
                    sublead_bjet_mass = choose_jet(dijets["second_jet"].mass, 0, -999.0)
                    sublead_bjet_charge = choose_jet(dijets["second_jet"].charge, 0, -999.0)
                    sublead_bjet_btagPNetB = choose_jet(dijets["second_jet"].btagPNetB, 0, -999.0)
                    sublead_bjet_PNetRegPtRawCorr = choose_jet(dijets["second_jet"].PNetRegPtRawCorr, 0, -999.0)
                    sublead_bjet_PNetRegPtRawCorrNeutrino = choose_jet(dijets["second_jet"].PNetRegPtRawCorrNeutrino, 0, -999.0)
                    sublead_bjet_PNetRegPtRawRes = choose_jet(dijets["second_jet"].PNetRegPtRawRes, 0, -999.0)

                    dijet_pt = choose_jet(dijets.pt, 0, -999.0)
                    dijet_eta = choose_jet(dijets.eta, 0, -999.0)
                    dijet_phi = choose_jet(dijets.phi, 0, -999.0)
                    dijet_mass = choose_jet(dijets.mass, 0, -999.0)
                    dijet_charge = choose_jet(dijets.charge, 0, -999.0)

                    # Get the HHbbgg object
                    HHbbgg = get_HHbbgg(self, diphotons, dijets)

                    # Write the variables in diphotons
                    diphotons["HHbbggCandidate_pt"] = HHbbgg.obj_HHbbgg.pt
                    diphotons["HHbbggCandidate_eta"] = HHbbgg.obj_HHbbgg.eta
                    diphotons["HHbbggCandidate_phi"] = HHbbgg.obj_HHbbgg.phi
                    diphotons["HHbbggCandidate_mass"] = HHbbgg.obj_HHbbgg.mass

                    diphotons["lead_bjet_pt"] = lead_bjet_pt
                    diphotons["lead_bjet_eta"] = lead_bjet_eta
                    diphotons["lead_bjet_phi"] = lead_bjet_phi
                    diphotons["lead_bjet_mass"] = lead_bjet_mass
                    diphotons["lead_bjet_charge"] = lead_bjet_charge
                    diphotons["lead_bjet_btagPNetB"] = lead_bjet_btagPNetB
                    diphotons["lead_bjet_PNetRegPtRawCorr"] = lead_bjet_PNetRegPtRawCorr
                    diphotons["lead_bjet_PNetRegPtRawCorrNeutrino"] = lead_bjet_PNetRegPtRawCorrNeutrino
                    diphotons["lead_bjet_PNetRegPtRawRes"] = lead_bjet_PNetRegPtRawRes

                    diphotons["sublead_bjet_pt"] = sublead_bjet_pt
                    diphotons["sublead_bjet_eta"] = sublead_bjet_eta
                    diphotons["sublead_bjet_phi"] = sublead_bjet_phi
                    diphotons["sublead_bjet_mass"] = sublead_bjet_mass
                    diphotons["sublead_bjet_charge"] = sublead_bjet_charge
                    diphotons["sublead_bjet_btagPNetB"] = sublead_bjet_btagPNetB
                    diphotons["sublead_bjet_PNetRegPtRawCorr"] = sublead_bjet_PNetRegPtRawCorr
                    diphotons["sublead_bjet_PNetRegPtRawCorrNeutrino"] = sublead_bjet_PNetRegPtRawCorrNeutrino
                    diphotons["sublead_bjet_PNetRegPtRawRes"] = sublead_bjet_PNetRegPtRawRes

                    diphotons["dijet_pt"] = dijet_pt
                    diphotons["dijet_eta"] = dijet_eta
                    diphotons["dijet_phi"] = dijet_phi
                    diphotons["dijet_mass"] = dijet_mass
                    diphotons["dijet_charge"] = dijet_charge

                    diphotons["pholead_PtOverM"] = HHbbgg.pho_lead.pt / HHbbgg.obj_diphoton.mass
                    diphotons["phosublead_PtOverM"] = HHbbgg.pho_sublead.pt / HHbbgg.obj_diphoton.mass

                    diphotons["FirstJet_PtOverM"] = diphotons["lead_bjet_pt"] / diphotons["dijet_mass"]
                    diphotons["SecondJet_PtOverM"] = diphotons["sublead_bjet_pt"] / diphotons["dijet_mass"]

                    diphotons["CosThetaStar_CS"] = getCosThetaStar_CS(HHbbgg, 6800)
                    diphotons["CosThetaStar_gg"] = getCosThetaStar_gg(HHbbgg)
                    diphotons["CosThetaStar_jj"] = getCosThetaStar_jj(HHbbgg)

                    diphotons["DeltaR_j1g1"] = DeltaR(HHbbgg.first_jet, HHbbgg.pho_lead)
                    diphotons["DeltaR_j2g1"] = DeltaR(HHbbgg.second_jet, HHbbgg.pho_lead)
                    diphotons["DeltaR_j1g2"] = DeltaR(HHbbgg.first_jet, HHbbgg.pho_sublead)
                    diphotons["DeltaR_j2g2"] = DeltaR(HHbbgg.second_jet, HHbbgg.pho_sublead)

                    DeltaR_comb = awkward.Array([diphotons["DeltaR_j1g1"], diphotons["DeltaR_j2g1"], diphotons["DeltaR_j1g2"], diphotons["DeltaR_j2g2"]])

                    diphotons["DeltaR_jg_min"] = awkward.min(DeltaR_comb, axis=0)

                    diphotons = diphotons[
                        diphotons["sublead_bjet_pt"] > -998
                    ]

                    # Addition of lepton info-> Taken from the top workflow. This part of the code was orignally written by Florain Mausolf
                    # Adding a 'generation' field to electrons and muons
                    sel_electrons['generation'] = awkward.ones_like(sel_electrons.pt)
                    sel_muons['generation'] = 2 * awkward.ones_like(sel_muons.pt)

                    # Combine electrons and muons into a single leptons collection
                    leptons = awkward.concatenate([sel_electrons, sel_muons], axis=1)
                    leptons = awkward.with_name(leptons, "PtEtaPhiMCandidate")

                    # Sort leptons by pt in descending order
                    leptons = leptons[awkward.argsort(leptons.pt, ascending=False)]

                    n_leptons = awkward.num(leptons)
                    diphotons["n_leptons"] = n_leptons

                    # Annotate diphotons with selected leptons properties
                    lepton_properties = ["pt", "eta", "phi", "mass", "charge", "generation"]
                    for i in range(self.num_leptons_to_store):  # Number of leptons to select
                        for prop in lepton_properties:
                            key = f"lepton{i+1}_{prop}"
                            # Retrieve the value using the choose_jet function (which can be used for leptons as well)
                            value = choose_jet(getattr(leptons, prop), i, -999.0)
                            # Store the value in the diphotons dictionary
                            diphotons[key] = value

                    # addition of fatjets and matching subjets, genjet
                    n_fatjets = awkward.num(fatjets)
                    diphotons["n_fatjets"] = n_fatjets

                    fatjet_properties = fatjets.fields
                    subjet_properties = subjets.fields
                    if self.data_kind == "mc":
                        genjetAK8_properties = genjetsAK8.fields

                    for i in range(self.num_fatjets_to_store):  # Number of fatjets to select
                        for prop in fatjet_properties:
                            if prop[-1] == "G":  # Few of the Idx variables are repeated with name ending with G (eg: 'subJetIdx1G'). Have to figure out why is this the case
                                continue
                            key = f"fatjet{i+1}_{prop}"
                            # Retrieve the value using the choose_jet function (which can be used for fatjets as well)
                            value = choose_jet(fatjets[prop], i, -999.0)

                            if prop == "genJetAK8Idx":  # add info of matched GenJetAK8
                                for prop_genJetAK8 in genjetAK8_properties:
                                    key_genJetAK8 = f"fatjet{i+1}_genjetAK8_{prop_genJetAK8}"
                                    # Retrieve the value using the choose_jet function (which can also be used here)
                                    value_genJetAK8 = choose_jet(genjetsAK8[prop_genJetAK8], value, -999.0)
                                    # Store the value in the diphotons dictionary
                                    diphotons[key_genJetAK8] = value_genJetAK8
                                continue  # not saving the index values

                            if prop in ["subJetIdx1", "subJetIdx2"]:  # add info of matched SubJets
                                subjet_name = prop.replace("Idx", "").lower()
                                for prop_subjet in subjet_properties:
                                    key_subjet = f"fatjet{i+1}_{subjet_name}_{prop_subjet}"
                                    # Retrieve the value using the choose_jet function (which can also be used here)
                                    value_subjet = choose_jet(subjets[prop_subjet], value, -999.0)
                                    # Store the value in the diphotons dictionary
                                    diphotons[key_subjet] = value_subjet
                                continue  # not saving the index values

                            # Store the value in the diphotons dictionary
                            diphotons[key] = value

                    # adding MET to parquet, adding all variables for now
                    puppiMET_properties = puppiMET.fields
                    for prop in puppiMET_properties:
                        key = f"puppiMET_{prop}"
                        # Retrieve the value using the choose_jet function (which can be used for puppiMET as well)
                        value = getattr(puppiMET, prop)
                        # Store the value in the diphotons dictionary
                        diphotons[key] = value

                    ## ----------  End of the HHTobbgg part ----------

                    # run taggers on the events list with added diphotons
                    # the shape here is ensured to be broadcastable
                    for tagger in self.taggers:
                        (
                            diphotons["_".join([tagger.name, str(tagger.priority)])],
                            tagger_extra,
                        ) = tagger(
                            events, diphotons
                        )  # creates new column in diphotons - tagger priority, or 0, also return list of histrograms here?
                        histos_etc.update(tagger_extra)

                    # if there are taggers to run, arbitrate by them first
                    # Deal with order of tagger priorities
                    # Turn from diphoton jagged array to whether or not an event was selected
                    if len(self.taggers):
                        counts = awkward.num(diphotons.pt, axis=1)
                        flat_tags = numpy.stack(
                            (
                                awkward.flatten(
                                    diphotons[
                                        "_".join([tagger.name, str(tagger.priority)])
                                    ]
                                )
                                for tagger in self.taggers
                            ),
                            axis=1,
                        )
                        tags = awkward.from_regular(
                            awkward.unflatten(flat_tags, counts), axis=2
                        )
                        winner = awkward.min(tags[tags != 0], axis=2)
                        diphotons["best_tag"] = winner

                        # lowest priority is most important (ascending sort)
                        # leave in order of diphoton pT in case of ties (stable sort)
                        sorted = awkward.argsort(diphotons.best_tag, stable=True)
                        diphotons = diphotons[sorted]

                    diphotons = awkward.firsts(diphotons)
                    # set diphotons as part of the event record
                    events[f"diphotons_{do_variation}"] = diphotons
                    # annotate diphotons with event information
                    diphotons["event"] = events.event
                    diphotons["lumi"] = events.luminosityBlock
                    diphotons["run"] = events.run
                    # nPV just for validation of pileup reweighting
                    diphotons["nPV"] = events.PV.npvs
                    diphotons["fixedGridRhoAll"] = events.Rho.fixedGridRhoAll
                    # annotate diphotons with dZ information (difference between z position of GenVtx and PV) as required by flashggfinalfits
                    if self.data_kind == "mc":
                        diphotons["genWeight"] = events.genWeight
                        diphotons["dZ"] = events.GenVtx.z - events.PV.z
                        # Necessary for differential xsec measurements in final fits ("truth" variables)
                        diphotons["HTXS_Higgs_pt"] = events.HTXS.Higgs_pt
                        diphotons["HTXS_Higgs_y"] = events.HTXS.Higgs_y
                        diphotons["HTXS_njets30"] = events.HTXS.njets30  # Need to clarify if this variable is suitable, does it fulfill abs(eta_j) < 2.5? Probably not
                        # Preparation for HTXS measurements later, start with stage 0 to disentangle VH into WH and ZH for final fits
                        diphotons["HTXS_stage_0"] = events.HTXS.stage_0
                    # Fill zeros for data because there is no GenVtx for data, obviously
                    else:
                        diphotons["dZ"] = awkward.zeros_like(events.PV.z)

                    # drop events without a preselected diphoton candidate
                    # drop events without a tag, if there are tags
                    if len(self.taggers):
                        selection_mask = ~(
                            awkward.is_none(diphotons)
                            | awkward.is_none(diphotons.best_tag)
                        )
                        diphotons = diphotons[selection_mask]
                    else:  # drop events without a preselected diphoton candidate
                        selection_mask = ~awkward.is_none(diphotons)
                        diphotons = diphotons[selection_mask]

                    # return if there is no surviving events
                    if len(diphotons) == 0:
                        logger.debug("No surviving events in this run, return now!")
                        return histos_etc
                    if self.data_kind == "mc":
                        # initiate Weight container here, after selection, since event selection cannot easily be applied to weight container afterwards
                        event_weights = Weights(size=len(events[selection_mask]))

                        # corrections to event weights:
                        for correction_name in correction_names:
                            if correction_name in available_weight_corrections:
                                logger.info(
                                    f"Adding correction {correction_name} to weight collection of dataset {dataset_name}"
                                )
                                varying_function = available_weight_corrections[
                                    correction_name
                                ]
                                event_weights = varying_function(
                                    events=events[selection_mask],
                                    photons=events[f"diphotons_{do_variation}"][
                                        selection_mask
                                    ],
                                    weights=event_weights,
                                    dataset_name=dataset_name,
                                    year=self.year[dataset_name][0],
                                )

                        # systematic variations of event weights go to nominal output dataframe:
                        if do_variation == "nominal":
                            for systematic_name in systematic_names:
                                if systematic_name in available_weight_systematics:
                                    logger.info(
                                        f"Adding systematic {systematic_name} to weight collection of dataset {dataset_name}"
                                    )
                                    if systematic_name == "LHEScale":
                                        if hasattr(events, "LHEScaleWeight"):
                                            diphotons["nweight_LHEScale"] = awkward.num(
                                                events.LHEScaleWeight[selection_mask],
                                                axis=1,
                                            )
                                            diphotons[
                                                "weight_LHEScale"
                                            ] = events.LHEScaleWeight[selection_mask]
                                        else:
                                            logger.info(
                                                f"No {systematic_name} Weights in dataset {dataset_name}"
                                            )
                                    elif systematic_name == "LHEPdf":
                                        if hasattr(events, "LHEPdfWeight"):
                                            # two AlphaS weights are removed
                                            diphotons["nweight_LHEPdf"] = (
                                                awkward.num(
                                                    events.LHEPdfWeight[selection_mask],
                                                    axis=1,
                                                )
                                                - 2
                                            )
                                            diphotons[
                                                "weight_LHEPdf"
                                            ] = events.LHEPdfWeight[selection_mask][
                                                :, :-2
                                            ]
                                        else:
                                            logger.info(
                                                f"No {systematic_name} Weights in dataset {dataset_name}"
                                            )
                                    else:
                                        varying_function = available_weight_systematics[
                                            systematic_name
                                        ]
                                        event_weights = varying_function(
                                            events=events[selection_mask],
                                            photons=events[f"diphotons_{do_variation}"][
                                                selection_mask
                                            ],
                                            weights=event_weights,
                                            dataset_name=dataset_name,
                                            year=self.year[dataset_name][0],
                                        )

                        diphotons["weight_central"] = event_weights.weight()
                        # Store variations with respect to central weight
                        if do_variation == "nominal":
                            if len(event_weights.variations):
                                logger.info(
                                    "Adding systematic weight variations to nominal output file."
                                )
                            for modifier in event_weights.variations:
                                diphotons["weight_" + modifier] = event_weights.weight(
                                    modifier=modifier
                                )

                        # Multiply weight by genWeight for normalisation in post-processing chain
                        event_weights._weight = (
                            events["genWeight"][selection_mask]
                            * diphotons["weight_central"]
                        )
                        diphotons["weight"] = event_weights.weight()

                    # Add weight variables (=1) for data for consistent datasets
                    else:
                        diphotons["weight_central"] = awkward.ones_like(
                            diphotons["event"]
                        )
                        diphotons["weight"] = awkward.ones_like(diphotons["event"])

                    ### Add mass resolution uncertainty
                    # Note that pt*cosh(eta) is equal to the energy of a four vector
                    # Note that you need to call it slightly different than in the output of HiggsDNA as pho_lead -> lead is only done in dumping utils

                    if (self.data_kind == "mc" and self.doFlow_corrections):
                        diphotons["sigma_m_over_m"] = 0.5 * numpy.sqrt(
                            (
                                diphotons["pho_lead"].raw_energyErr
                                / (
                                    diphotons["pho_lead"].pt
                                    * numpy.cosh(diphotons["pho_lead"].eta)
                                )
                            )
                            ** 2
                            + (
                                diphotons["pho_sublead"].raw_energyErr
                                / (
                                    diphotons["pho_sublead"].pt
                                    * numpy.cosh(diphotons["pho_sublead"].eta)
                                )
                            )
                            ** 2
                        )

                        diphotons["sigma_m_over_m_corr"] = 0.5 * numpy.sqrt(
                            (
                                diphotons["pho_lead"].energyErr
                                / (
                                    diphotons["pho_lead"].pt
                                    * numpy.cosh(diphotons["pho_lead"].eta)
                                )
                            )
                            ** 2
                            + (
                                diphotons["pho_sublead"].energyErr
                                / (
                                    diphotons["pho_sublead"].pt
                                    * numpy.cosh(diphotons["pho_sublead"].eta)
                                )
                            )
                            ** 2
                        )

                    else:
                        diphotons["sigma_m_over_m"] = 0.5 * numpy.sqrt(
                            (
                                diphotons["pho_lead"].energyErr
                                / (
                                    diphotons["pho_lead"].pt
                                    * numpy.cosh(diphotons["pho_lead"].eta)
                                )
                            )
                            ** 2
                            + (
                                diphotons["pho_sublead"].energyErr
                                / (
                                    diphotons["pho_sublead"].pt
                                    * numpy.cosh(diphotons["pho_sublead"].eta)
                                )
                            )
                            ** 2
                        )

                    # This is the mass SigmaM/M value including the smearing term from the Scale and smearing
                    # The implementation follows the flashGG implementation -> https://github.com/cms-analysis/flashgg/blob/4edea8897e2a4b0518dca76ba6c9909c20c40ae7/DataFormats/src/Photon.cc#L293
                    # adittional flashGG link when the smearing of the SigmaE/E smearing is called -> https://github.com/cms-analysis/flashgg/blob/4edea8897e2a4b0518dca76ba6c9909c20c40ae7/Systematics/plugins/PhotonSigEoverESmearingEGMTool.cc#L83C40-L83C45
                    # Just a reminder, the pt/energy of teh data is not smearing, but the smearing term is added to the data sigma_m_over_m
                    if (self.Smear_sigma_m):

                        if (self.doFlow_corrections and self.data_kind == "mc"):
                            # Adding the smeared BDT error to the ntuples!
                            diphotons["pho_lead","energyErr_Smeared"] = numpy.sqrt((diphotons["pho_lead"].raw_energyErr)**2 + (diphotons["pho_lead"].rho_smear * ((diphotons["pho_lead"].pt * numpy.cosh(diphotons["pho_lead"].eta)))) ** 2)
                            diphotons["pho_sublead","energyErr_Smeared"] = numpy.sqrt((diphotons["pho_sublead"].raw_energyErr) ** 2 + (diphotons["pho_sublead"].rho_smear * ((diphotons["pho_sublead"].pt * numpy.cosh(diphotons["pho_sublead"].eta)))) ** 2)

                            diphotons["sigma_m_over_m_Smeared"] = 0.5 * numpy.sqrt(
                                (
                                    numpy.sqrt((diphotons["pho_lead"].raw_energyErr)**2 + (diphotons["pho_lead"].rho_smear * ((diphotons["pho_lead"].pt * numpy.cosh(diphotons["pho_lead"].eta)))) ** 2)
                                    / (
                                        diphotons["pho_lead"].pt
                                        * numpy.cosh(diphotons["pho_lead"].eta)
                                    )
                                )
                                ** 2
                                + (
                                    numpy.sqrt((diphotons["pho_sublead"].raw_energyErr) ** 2 + (diphotons["pho_sublead"].rho_smear * ((diphotons["pho_sublead"].pt * numpy.cosh(diphotons["pho_sublead"].eta)))) ** 2)
                                    / (
                                        diphotons["pho_sublead"].pt
                                        * numpy.cosh(diphotons["pho_sublead"].eta)
                                    )
                                )
                                ** 2
                            )

                            diphotons["sigma_m_over_m_Smeared_Corr"] = 0.5 * numpy.sqrt(
                                (
                                    numpy.sqrt((diphotons["pho_lead"].energyErr)**2 + (diphotons["pho_lead"].rho_smear * ((diphotons["pho_lead"].pt * numpy.cosh(diphotons["pho_lead"].eta)))) ** 2)
                                    / (
                                        diphotons["pho_lead"].pt
                                        * numpy.cosh(diphotons["pho_lead"].eta)
                                    )
                                )
                                ** 2
                                + (
                                    numpy.sqrt((diphotons["pho_sublead"].energyErr) ** 2 + (diphotons["pho_sublead"].rho_smear * ((diphotons["pho_sublead"].pt * numpy.cosh(diphotons["pho_sublead"].eta)))) ** 2)
                                    / (
                                        diphotons["pho_sublead"].pt
                                        * numpy.cosh(diphotons["pho_sublead"].eta)
                                    )
                                )
                                ** 2
                            )

                        else:
                            # Adding the smeared BDT error to the ntuples!
                            diphotons["pho_lead","energyErr_Smeared"] = numpy.sqrt((diphotons["pho_lead"].energyErr)**2 + (diphotons["pho_lead"].rho_smear * ((diphotons["pho_lead"].pt * numpy.cosh(diphotons["pho_lead"].eta)))) ** 2)
                            diphotons["pho_sublead","energyErr_Smeared"] = numpy.sqrt((diphotons["pho_sublead"].energyErr) ** 2 + (diphotons["pho_sublead"].rho_smear * ((diphotons["pho_sublead"].pt * numpy.cosh(diphotons["pho_sublead"].eta)))) ** 2)

                            diphotons["sigma_m_over_m_Smeared"] = 0.5 * numpy.sqrt(
                                (
                                    numpy.sqrt((diphotons["pho_lead"].energyErr)**2 + (diphotons["pho_lead"].rho_smear * ((diphotons["pho_lead"].pt * numpy.cosh(diphotons["pho_lead"].eta)))) ** 2)
                                    / (
                                        diphotons["pho_lead"].pt
                                        * numpy.cosh(diphotons["pho_lead"].eta)
                                    )
                                )
                                ** 2
                                + (
                                    numpy.sqrt((diphotons["pho_sublead"].energyErr) ** 2 + (diphotons["pho_sublead"].rho_smear * ((diphotons["pho_sublead"].pt * numpy.cosh(diphotons["pho_sublead"].eta)))) ** 2)
                                    / (
                                        diphotons["pho_sublead"].pt
                                        * numpy.cosh(diphotons["pho_sublead"].eta)
                                    )
                                )
                                ** 2
                            )

                    # Decorrelating the mass resolution - Still need to supress the decorrelator noises
                    if self.doDeco:

                        # Decorrelate nominal sigma_m_over_m
                        diphotons["sigma_m_over_m_nominal_decorr"] = decorrelate_mass_resolution(diphotons, type="nominal", year=self.year[dataset_name][0])

                        # decorrelate smeared nominal sigma_m_overm_m
                        if (self.Smear_sigma_m):
                            diphotons["sigma_m_over_m_smeared_decorr"] = decorrelate_mass_resolution(diphotons, type="smeared", year=self.year[dataset_name][0])

                        # decorrelate flow corrected sigma_m_over_m
                        if (self.doFlow_corrections):
                            diphotons["sigma_m_over_m_corr_decorr"] = decorrelate_mass_resolution(diphotons, type="corr", year=self.year[dataset_name][0])

                        # decorrelate flow corrected smeared sigma_m_over_m
                        if (self.doFlow_corrections and self.Smear_sigma_m):
                            diphotons["sigma_m_over_m_corr_smeared_decorr"] = decorrelate_mass_resolution(diphotons, type="corr_smeared", year=self.year[dataset_name][0])

                        # Instead of the nominal sigma_m_over_m, we will use the smeared version of it -> (https://indico.cern.ch/event/1319585/#169-update-on-the-run-3-mass-r)
                        # else:
                        #    warnings.warn("Smeamering need to be applied in order to decorrelate the (Smeared) mass resolution. -- Exiting!")
                        #    sys.exit(0)

                    if self.output_location is not None:
                        if self.output_format == "root":
                            df = diphoton_list_to_pandas(self, diphotons)
                        else:
                            akarr = diphoton_ak_array(self, diphotons)

                            # Remove fixedGridRhoAll from photons to avoid having event-level info per photon
                            akarr = akarr[
                                [
                                    field
                                    for field in akarr.fields
                                    if "lead_fixedGridRhoAll" not in field
                                ]
                            ]

                        fname = (
                            events.behavior[
                                "__events_factory__"
                            ]._partition_key.replace("/", "_")
                            + ".%s" % self.output_format
                        )
                        subdirs = []
                        if "dataset" in events.metadata:
                            subdirs.append(events.metadata["dataset"])
                        subdirs.append(do_variation)
                        if self.output_format == "root":
                            dump_pandas(self, df, fname, self.output_location, subdirs)
                        else:
                            dump_ak_array(
                                self, akarr, fname, self.output_location, metadata, subdirs,
                            )

        return histos_etc

    def postprocess(self, accumulant: Dict[Any, Any]) -> Any:
        pass
