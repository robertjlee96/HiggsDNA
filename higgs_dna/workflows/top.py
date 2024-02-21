from higgs_dna.workflows.base import HggBaseProcessor

from higgs_dna.tools.SC_eta import add_photon_SC_eta
from higgs_dna.tools.EELeak_region import veto_EEleak_flag
from higgs_dna.selections.photon_selections import photon_preselection
from higgs_dna.selections.lepton_selections import select_electrons, select_muons
from higgs_dna.selections.jet_selections import select_jets, jetvetomap
from higgs_dna.selections.lumi_selections import select_lumis
from higgs_dna.utils.dumping_utils import diphoton_ak_array, dump_ak_array, diphoton_list_to_pandas, dump_pandas
from higgs_dna.utils.misc_utils import choose_jet
from higgs_dna.tools.flow_corrections import calculate_flow_corrections

from higgs_dna.systematics import object_systematics as available_object_systematics
from higgs_dna.systematics import object_corrections as available_object_corrections
from higgs_dna.systematics import weight_systematics as available_weight_systematics
from higgs_dna.systematics import weight_corrections as available_weight_corrections

import functools
import warnings
from typing import Any, Dict, List, Optional
import awkward as ak
import numpy
import vector
from coffea.analysis_tools import Weights
from copy import deepcopy

import logging

logger = logging.getLogger(__name__)

vector.register_awkward()


class TopProcessor(HggBaseProcessor):  # type: ignore

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

    def process_extra(self, events: ak.Array) -> ak.Array:
        return events, {}

    def postprocess(self, accumulant: Dict[Any, Any]) -> Any:
        pass

    def process(self, events: ak.Array) -> Dict[Any, Any]:

        print("\n \t INFO: running top processor. \n")
        dataset_name = events.metadata["dataset"]

        # data or monte carlo?
        self.data_kind = "mc" if hasattr(events, "GenPart") else "data"

        # here we start recording possible coffea accumulators
        # most likely histograms, could be counters, arrays, ...
        histos_etc = {}
        histos_etc[dataset_name] = {}
        if self.data_kind == "mc":
            histos_etc[dataset_name]["nTot"] = int(
                ak.num(events.genWeight, axis=0)
            )
            histos_etc[dataset_name]["nPos"] = int(ak.sum(events.genWeight > 0))
            histos_etc[dataset_name]["nNeg"] = int(ak.sum(events.genWeight < 0))
            histos_etc[dataset_name]["nEff"] = int(
                histos_etc[dataset_name]["nPos"] - histos_etc[dataset_name]["nNeg"]
            )
            histos_etc[dataset_name]["genWeightSum"] = float(
                ak.sum(events.genWeight)
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
            metadata["sum_genw_presel"] = str(ak.sum(events.genWeight))
        else:
            metadata["sum_genw_presel"] = "Data"

        # apply filters and triggers
        events = self.apply_filters_and_triggers(events)

        # we need ScEta for corrections and systematics, which is not present in NanoAODv11 but can be calculated using PV
        events.Photon = add_photon_SC_eta(events.Photon, events.PV)

        # add veto EE leak branch for photons, could also be used for electrons
        if self.year[dataset_name][0] == "2022postEE":
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

                    # Computing the normalizing flow correction
                    if self.data_kind == "mc" and self.doFlow_corrections:

                        # Applying the Flow corrections to all photons before pre-selection
                        counts = ak.num(photons)
                        corrected_inputs,var_list = calculate_flow_corrections(photons, events, self.meta["flashggPhotons"]["flow_inputs"], self.meta["flashggPhotons"]["Isolation_transform_order"], year=self.year[dataset_name][0])

                        # Store the raw nanoAOD value and update photon ID MVA value for preselection
                        photons["mvaID_run3"] = ak.unflatten(self.add_photonid_mva_run3(photons, events), counts)
                        photons["mvaID_nano"] = photons["mvaID"]

                        # Store the raw values of the inputs and update the input values with the corrections since some variables used in the preselection
                        for i in range(len(var_list)):
                            photons["raw_" + str(var_list[i])] = photons[str(var_list[i])]
                            photons[str(var_list[i])] = ak.unflatten(corrected_inputs[:,i] , counts)

                        photons["mvaID"] = ak.unflatten(self.add_photonid_mva_run3(photons, events), counts)

                    # photon preselection
                    photons = photon_preselection(self, photons, events, year=self.year[dataset_name][0])
                    # sort photons in each event descending in pt
                    # make descending-pt combinations of photons
                    photons = photons[ak.argsort(photons.pt, ascending=False)]
                    photons["charge"] = ak.zeros_like(
                        photons.pt
                    )  # added this because charge is not a property of photons in nanoAOD v11. We just assume every photon has charge zero...
                    diphotons = ak.combinations(
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
                    diphotons = ak.with_name(diphotons, "PtEtaPhiMCandidate")

                    # sort diphotons by pT
                    diphotons = diphotons[
                        ak.argsort(diphotons.pt, ascending=False)
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
                    jets = ak.zip(
                        {
                            "pt": Jets.pt,
                            "eta": Jets.eta,
                            "phi": Jets.phi,
                            "mass": Jets.mass,
                            "charge": ak.zeros_like(
                                Jets.pt
                            ),  # added this because jet charge is not a property of photons in nanoAOD v11. We just need the charge to build jet collection.
                            "hFlav": Jets.hadronFlavour
                            if self.data_kind == "mc"
                            else ak.zeros_like(Jets.pt),
                            "btagDeepFlav_B": Jets.btagDeepFlavB,
                            "btagDeepFlav_CvB": Jets.btagDeepFlavCvB,
                            "btagDeepFlav_CvL": Jets.btagDeepFlavCvL,
                            "btagDeepFlav_QG": Jets.btagDeepFlavQG,
                        }
                    )
                    jets = ak.with_name(jets, "PtEtaPhiMCandidate")

                    electrons = ak.zip(
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
                    electrons = ak.with_name(electrons, "PtEtaPhiMCandidate")

                    muons = ak.zip(
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
                    muons = ak.with_name(muons, "PtEtaPhiMCandidate")

                    # lepton cleaning
                    sel_electrons = electrons[
                        select_electrons(self, electrons, diphotons)
                    ]
                    sel_muons = muons[select_muons(self, muons, diphotons)]

                    # jet selection and pt ordering
                    jets = jets[
                        select_jets(self, jets, diphotons, sel_muons, sel_electrons)
                    ]
                    jets = jets[ak.argsort(jets.pt, ascending=False)]

                    # adding selected jets to events to be used in ctagging SF calculation
                    events["sel_jets"] = jets
                    n_jets = ak.num(jets)

                    num_jets = 6
                    jet_properties = ["pt", "eta", "phi", "mass", "charge", "btagDeepFlav_B"]
                    for i in range(num_jets):
                        for prop in jet_properties:
                            key = f"jet{i+1}_{prop}"
                            value = choose_jet(getattr(jets, prop), i, -999.0)
                            # Store the value in the diphotons dictionary
                            diphotons[key] = value
                    diphotons["n_jets"] = n_jets

                    # Adding a 'generation' field to electrons and muons
                    sel_electrons['generation'] = ak.ones_like(sel_electrons.pt)
                    sel_muons['generation'] = 2 * ak.ones_like(sel_muons.pt)

                    # Combine electrons and muons into a single leptons collection
                    leptons = ak.concatenate([sel_electrons, sel_muons], axis=1)
                    leptons = ak.with_name(leptons, "PtEtaPhiMCandidate")

                    # Sort leptons by pt in descending order
                    leptons = leptons[ak.argsort(leptons.pt, ascending=False)]

                    n_leptons = ak.num(leptons)
                    diphotons["n_leptons"] = n_leptons

                    # Annotate diphotons with selected leptons properties
                    lepton_properties = ["pt", "eta", "phi", "mass", "charge", "generation"]
                    num_leptons = 2  # Number of leptons to select
                    for i in range(num_leptons):
                        for prop in lepton_properties:
                            key = f"lepton{i+1}_{prop}"
                            # Retrieve the value using the choose_jet function (which can be used for leptons as well)
                            value = choose_jet(getattr(leptons, prop), i, -999.0)
                            # Store the value in the diphotons dictionary
                            diphotons[key] = value

                    diphotons = ak.firsts(diphotons)
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
                        diphotons["HTXS_Higgs_pt"] = events.HTXS.Higgs_pt
                        diphotons["HTXS_Higgs_y"] = events.HTXS.Higgs_y
                        diphotons["HTXS_njets30"] = events.HTXS.njets30
                        diphotons["HTXS_stage_0"] = events.HTXS.stage_0
                    else:
                        diphotons["dZ"] = ak.zeros_like(events.PV.z)

                    # drop events without a preselected diphoton candidate
                    selection_mask = ~ak.is_none(diphotons)
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
                                            diphotons["nweight_LHEScale"] = ak.num(
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
                                                ak.num(
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

                        if ak.num(events.LHEReweightingWeight)[0] > 0:
                            diphotons["LHEReweightingWeight"] = events.LHEReweightingWeight[selection_mask]
                            diphotons["LHEWeight"] = events.LHEWeight[selection_mask]

                    # Add weight variables (=1) for data for consistent datasets
                    else:
                        diphotons["weight_central"] = ak.ones_like(
                            diphotons["event"]
                        )
                        diphotons["weight"] = ak.ones_like(diphotons["event"])

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
