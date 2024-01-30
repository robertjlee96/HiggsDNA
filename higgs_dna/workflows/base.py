from higgs_dna.tools.chained_quantile import ChainedQuantileRegression
from higgs_dna.tools.diphoton_mva import calculate_diphoton_mva
from higgs_dna.tools.xgb_loader import load_bdt
from higgs_dna.tools.photonid_mva import calculate_photonid_mva, load_photonid_mva
from higgs_dna.tools.photonid_mva import calculate_photonid_mva_run3, load_photonid_mva_run3
from higgs_dna.tools.SC_eta import add_photon_SC_eta
from higgs_dna.tools.EELeak_region import veto_EEleak_flag
from higgs_dna.selections.photon_selections import photon_preselection
from higgs_dna.selections.lepton_selections import select_electrons, select_muons
from higgs_dna.selections.jet_selections import select_jets, jetvetomap
from higgs_dna.selections.lumi_selections import select_lumis
from higgs_dna.utils.dumping_utils import diphoton_ak_array, dump_ak_array, diphoton_list_to_pandas, dump_pandas
from higgs_dna.utils.misc_utils import choose_jet
from higgs_dna.tools.flow_corrections import calculate_flow_corrections

from higgs_dna.tools.mass_decorrelator import decorrelate_mass_resolution

# from higgs_dna.utils.dumping_utils import diphoton_list_to_pandas, dump_pandas
from higgs_dna.metaconditions import photon_id_mva_weights
from higgs_dna.metaconditions import diphoton as diphoton_mva_dir
from higgs_dna.systematics import object_systematics as available_object_systematics
from higgs_dna.systematics import object_corrections as available_object_corrections
from higgs_dna.systematics import weight_systematics as available_weight_systematics
from higgs_dna.systematics import weight_corrections as available_weight_corrections

import functools
import operator
import os
import warnings
from typing import Any, Dict, List, Optional
import awkward
import numpy
import sys
import vector
from coffea import processor
from coffea.analysis_tools import Weights
from copy import deepcopy

import logging

logger = logging.getLogger(__name__)

vector.register_awkward()


class HggBaseProcessor(processor.ProcessorABC):  # type: ignore
    def __init__(
        self,
        metaconditions: Dict[str, Any],
        systematics: Optional[Dict[str, List[str]]],
        corrections: Optional[Dict[str, List[str]]],
        apply_trigger: bool,
        output_location: Optional[str],
        taggers: Optional[List[Any]],
        trigger_group: str,
        analysis: str,
        skipCQR: bool,
        skipJetVetoMap: bool,
        year: Optional[Dict[str, List[str]]],
        fiducialCuts: str,
        doDeco: bool,
        Smear_sigma_m: bool,
        doFlow_corrections: bool,
        output_format: str,
    ) -> None:
        self.meta = metaconditions
        self.systematics = systematics if systematics is not None else {}
        self.corrections = corrections if corrections is not None else {}
        self.apply_trigger = apply_trigger
        self.output_location = output_location
        self.trigger_group = trigger_group
        self.analysis = analysis
        self.skipCQR = skipCQR
        self.skipJetVetoMap = skipJetVetoMap
        self.year = year if year is not None else {}
        self.fiducialCuts = fiducialCuts
        self.doDeco = doDeco
        self.Smear_sigma_m = Smear_sigma_m
        self.doFlow_corrections = doFlow_corrections
        self.output_format = output_format

        # muon selection cuts
        self.muon_pt_threshold = 10
        self.muon_max_eta = 2.4
        self.mu_iso_wp = "medium"
        self.global_muon = False

        # electron selection cuts
        self.electron_pt_threshold = 15
        self.electron_max_eta = 2.5
        self.el_iso_wp = "WP80"

        # jet selection cuts
        self.jet_dipho_min_dr = 0.4
        self.jet_pho_min_dr = 0.4
        self.jet_ele_min_dr = 0.4
        self.jet_muo_min_dr = 0.4
        self.jet_pt_threshold = 20
        self.jet_max_eta = 4.7

        self.clean_jet_dipho = True
        self.clean_jet_pho = True
        self.clean_jet_ele = False
        self.clean_jet_muo = False

        # diphoton preselection cuts
        self.min_pt_photon = 25.0
        self.min_pt_lead_photon = 35.0
        self.min_mvaid = -0.9
        self.max_sc_eta = 2.5
        self.gap_barrel_eta = 1.4442
        self.gap_endcap_eta = 1.566
        self.max_hovere = 0.08
        self.min_full5x5_r9 = 0.8
        self.max_chad_iso = 20.0
        self.max_chad_rel_iso = 0.3

        self.min_full5x5_r9_EB_high_r9 = 0.85
        self.min_full5x5_r9_EE_high_r9 = 0.9
        self.min_full5x5_r9_EB_low_r9 = 0.5
        self.min_full5x5_r9_EE_low_r9 = 0.8
        self.max_trkSumPtHollowConeDR03_EB_low_r9 = (
            6.0  # for v11, we cut on Photon_pfChargedIsoPFPV
        )
        self.max_trkSumPtHollowConeDR03_EE_low_r9 = 6.0  # Leaving the names of the preselection cut variables the same to change as little as possible
        self.max_sieie_EB_low_r9 = 0.015
        self.max_sieie_EE_low_r9 = 0.035
        self.max_pho_iso_EB_low_r9 = 4.0
        self.max_pho_iso_EE_low_r9 = 4.0

        self.eta_rho_corr = 1.5
        self.low_eta_rho_corr = 0.16544
        self.high_eta_rho_corr = 0.13212
        # EA values for Run3 from Egamma
        self.EA1_EB1 = 0.102056
        self.EA2_EB1 = -0.000398112
        self.EA1_EB2 = 0.0820317
        self.EA2_EB2 = -0.000286224
        self.EA1_EE1 = 0.0564915
        self.EA2_EE1 = -0.000248591
        self.EA1_EE2 = 0.0428606
        self.EA2_EE2 = -0.000171541
        self.EA1_EE3 = 0.0395282
        self.EA2_EE3 = -0.000121398
        self.EA1_EE4 = 0.0369761
        self.EA2_EE4 = -8.10369e-05
        self.EA1_EE5 = 0.0369417
        self.EA2_EE5 = -2.76885e-05
        self.e_veto = 0.5

        logger.debug(f"Setting up processor with metaconditions: {self.meta}")

        self.taggers = []
        if taggers is not None:
            self.taggers = taggers
            self.taggers.sort(key=lambda x: x.priority)

        self.prefixes = {"pho_lead": "lead", "pho_sublead": "sublead"}

        if not self.doDeco:
            logger.info("Skipping Mass resolution decorrelation as required")
        else:
            logger.info("Performing Mass resolution decorrelation as required")

        # build the chained quantile regressions
        if not self.skipCQR:
            try:
                self.chained_quantile: Optional[
                    ChainedQuantileRegression
                ] = ChainedQuantileRegression(**self.meta["PhoIdInputCorrections"])
            except Exception as e:
                warnings.warn(f"Could not instantiate ChainedQuantileRegression: {e}")
                self.chained_quantile = None
        else:
            logger.info("Skipping CQR as required")
            self.chained_quantile = None

        # initialize photonid_mva
        photon_id_mva_dir = os.path.dirname(photon_id_mva_weights.__file__)
        try:
            logger.debug(
                f"Looking for {self.meta['flashggPhotons']['photonIdMVAweightfile_EB']} in {photon_id_mva_dir}"
            )
            self.photonid_mva_EB = load_photonid_mva(
                os.path.join(
                    photon_id_mva_dir,
                    self.meta["flashggPhotons"]["photonIdMVAweightfile_EB"],
                )
            )
            self.photonid_mva_EE = load_photonid_mva(
                os.path.join(
                    photon_id_mva_dir,
                    self.meta["flashggPhotons"]["photonIdMVAweightfile_EE"],
                )
            )
        except Exception as e:
            warnings.warn(f"Could not instantiate PhotonID MVA on the fly: {e}")
            self.photonid_mva_EB = None
            self.photonid_mva_EE = None

        # initialize diphoton mva
        diphoton_weights_dir = os.path.dirname(diphoton_mva_dir.__file__)
        logger.debug(
            f"Base path to look for IDMVA weight files: {diphoton_weights_dir}"
        )

        try:
            self.diphoton_mva = load_bdt(
                os.path.join(
                    diphoton_weights_dir, self.meta["flashggDiPhotonMVA"]["weightFile"]
                )
            )
        except Exception as e:
            warnings.warn(f"Could not instantiate diphoton MVA: {e}")
            self.diphoton_mva = None

    def process_extra(self, events: awkward.Array) -> awkward.Array:
        raise NotImplementedError

    def apply_filters_and_triggers(self, events: awkward.Array) -> awkward.Array:
        # met filters
        met_filters = self.meta["flashggMetFilters"][self.data_kind]
        filtered = functools.reduce(
            operator.and_,
            (events.Flag[metfilter.split("_")[-1]] for metfilter in met_filters),
        )

        triggered = awkward.ones_like(filtered)
        if self.apply_trigger:
            trigger_names = []
            triggers = self.meta["TriggerPaths"][self.trigger_group][self.analysis]
            hlt = events.HLT
            for trigger in triggers:
                actual_trigger = trigger.replace("HLT_", "").replace("*", "")
                for field in hlt.fields:
                    if field.startswith(actual_trigger):
                        trigger_names.append(field)
            triggered = functools.reduce(
                operator.or_, (hlt[trigger_name] for trigger_name in trigger_names)
            )

        return events[filtered & triggered]

    def process(self, events: awkward.Array) -> Dict[Any, Any]:
        dataset_name = events.metadata["dataset"]

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

                    # Computing the normalizinf flow correction
                    if self.data_kind == "mc" and self.doFlow_corrections:

                        # Applyting the Flow corrections to all photons before pre-selection
                        counts = awkward.num(photons)
                        corrected_inputs,var_list = calculate_flow_corrections(photons, events, self.meta["flashggPhotons"]["flow_inputs"], self.meta["flashggPhotons"]["Isolation_transform_order"], year=self.year[dataset_name][0])

                        # adding the corrected values to the photons
                        for i in range(len(var_list)):
                            photons["corr_" + str(var_list[i])] = awkward.unflatten(corrected_inputs[:,i] , counts)

                        # Now adding the corrected mvaID to the photon entries
                        photons["mvaID_run3"] = awkward.unflatten(self.add_photonid_mva_run3(photons, events), counts)
                        photons["corr_mvaID_run3"] = awkward.unflatten(self.add_corr_photonid_mva_run3(photons, events), counts)

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

                    # adding selected jets to events to be used in ctagging SF calculation
                    events["sel_jets"] = jets
                    n_jets = awkward.num(jets)

                    first_jet_pt = choose_jet(jets.pt, 0, -999.0)
                    first_jet_eta = choose_jet(jets.eta, 0, -999.0)
                    first_jet_phi = choose_jet(jets.phi, 0, -999.0)
                    first_jet_mass = choose_jet(jets.mass, 0, -999.0)
                    first_jet_charge = choose_jet(jets.charge, 0, -999.0)

                    second_jet_pt = choose_jet(jets.pt, 1, -999.0)
                    second_jet_eta = choose_jet(jets.eta, 1, -999.0)
                    second_jet_phi = choose_jet(jets.phi, 1, -999.0)
                    second_jet_mass = choose_jet(jets.mass, 1, -999.0)
                    second_jet_charge = choose_jet(jets.charge, 1, -999.0)

                    diphotons["first_jet_pt"] = first_jet_pt
                    diphotons["first_jet_eta"] = first_jet_eta
                    diphotons["first_jet_phi"] = first_jet_phi
                    diphotons["first_jet_mass"] = first_jet_mass
                    diphotons["first_jet_charge"] = first_jet_charge

                    diphotons["second_jet_pt"] = second_jet_pt
                    diphotons["second_jet_eta"] = second_jet_eta
                    diphotons["second_jet_phi"] = second_jet_phi
                    diphotons["second_jet_mass"] = second_jet_mass
                    diphotons["second_jet_charge"] = second_jet_charge

                    diphotons["n_jets"] = n_jets

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
                    else:
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

                    if (self.data_kind == "mc" and self.doFlow_corrections):

                        diphotons["sigma_m_over_m_corr"] = 0.5 * numpy.sqrt(
                            (
                                diphotons["pho_lead"].corr_energyErr
                                / (
                                    diphotons["pho_lead"].pt
                                    * numpy.cosh(diphotons["pho_lead"].eta)
                                )
                            )
                            ** 2
                            + (
                                diphotons["pho_sublead"].corr_energyErr
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

                        if (self.data_kind == "mc" and self.doFlow_corrections):

                            diphotons["sigma_m_over_m_Smeared_corr"] = 0.5 * numpy.sqrt(
                                (
                                    numpy.sqrt((diphotons["pho_lead"].corr_energyErr)**2 + (diphotons["pho_lead"].rho_smear * ((diphotons["pho_lead"].pt * numpy.cosh(diphotons["pho_lead"].eta)))) ** 2)
                                    / (
                                        diphotons["pho_lead"].pt
                                        * numpy.cosh(diphotons["pho_lead"].eta)
                                    )
                                )
                                ** 2
                                + (
                                    numpy.sqrt((diphotons["pho_sublead"].corr_energyErr) ** 2 + (diphotons["pho_sublead"].rho_smear * ((diphotons["pho_sublead"].pt * numpy.cosh(diphotons["pho_sublead"].eta)))) ** 2)
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
        raise NotImplementedError

    def add_diphoton_mva(
        self, diphotons: awkward.Array, events: awkward.Array
    ) -> awkward.Array:
        return calculate_diphoton_mva(
            (self.diphoton_mva, self.meta["flashggDiPhotonMVA"]["inputs"]),
            diphotons,
            events,
        )

    def add_photonid_mva(
        self, photons: awkward.Array, events: awkward.Array
    ) -> awkward.Array:
        photons["fixedGridRhoAll"] = events.Rho.fixedGridRhoAll * awkward.ones_like(
            photons.pt
        )
        counts = awkward.num(photons, axis=-1)
        photons = awkward.flatten(photons)
        isEB = awkward.to_numpy(numpy.abs(photons.eta) < 1.5)
        mva_EB = calculate_photonid_mva(
            (self.photonid_mva_EB, self.meta["flashggPhotons"]["inputs_EB"]), photons
        )
        mva_EE = calculate_photonid_mva(
            (self.photonid_mva_EE, self.meta["flashggPhotons"]["inputs_EE"]), photons
        )
        mva = awkward.where(isEB, mva_EB, mva_EE)
        photons["mvaID"] = mva

        return awkward.unflatten(photons, counts)

    def add_photonid_mva_run3(
        self, photons: awkward.Array, events: awkward.Array
    ) -> awkward.Array:

        preliminary_path = os.path.join(os.path.dirname(__file__), '../tools/flows/run3_mvaID_models/')
        photonid_mva_EB, photonid_mva_EE = load_photonid_mva_run3(preliminary_path)

        rho = events.Rho.fixedGridRhoAll * awkward.ones_like(photons.pt)
        rho = awkward.flatten(rho)

        photons = awkward.flatten(photons)

        isEB = awkward.to_numpy(numpy.abs(photons.eta) < 1.5)
        mva_EB = calculate_photonid_mva_run3(
            [photonid_mva_EB, self.meta["flashggPhotons"]["inputs_EB"]], photons , rho
        )
        mva_EE = calculate_photonid_mva_run3(
            [photonid_mva_EE, self.meta["flashggPhotons"]["inputs_EE"]], photons, rho
        )
        mva = awkward.where(isEB, mva_EB, mva_EE)
        photons["mvaID_run3"] = mva

        return mva

    def add_corr_photonid_mva_run3(
        self, photons: awkward.Array, events: awkward.Array
    ) -> awkward.Array:

        preliminary_path = os.path.join(os.path.dirname(__file__), '../tools/flows/run3_mvaID_models/')
        photonid_mva_EB, photonid_mva_EE = load_photonid_mva_run3(preliminary_path)

        rho = events.Rho.fixedGridRhoAll * awkward.ones_like(photons.pt)
        rho = awkward.flatten(rho)

        photons = awkward.flatten(photons)

        # Now calculating the corrected mvaID
        isEB = awkward.to_numpy(numpy.abs(photons.eta) < 1.5)
        corr_mva_EB = calculate_photonid_mva_run3(
            [photonid_mva_EB, self.meta["flashggPhotons"]["inputs_EB_corr"]], photons, rho
        )
        corr_mva_EE = calculate_photonid_mva_run3(
            [photonid_mva_EE, self.meta["flashggPhotons"]["inputs_EE_corr"]], photons, rho
        )
        corr_mva = awkward.where(isEB, corr_mva_EB, corr_mva_EE)

        return corr_mva
