from higgs_dna.workflows.base import HggBaseProcessor
from higgs_dna.systematics import object_systematics as available_object_systematics
from higgs_dna.systematics import object_corrections as available_object_corrections
from higgs_dna.systematics import weight_systematics as available_weight_systematics
from higgs_dna.systematics import weight_corrections as available_weight_corrections
from higgs_dna.selections.photon_selections import photon_preselection
from higgs_dna.selections.lumi_selections import select_lumis
from higgs_dna.utils.dumping_utils import diphoton_list_to_pandas, dump_pandas
from higgs_dna.tools.SC_eta import add_photon_SC_eta
from higgs_dna.tools.flow_corrections import calculate_flow_corrections
from typing import Any, Dict, List, Optional
import awkward as ak
import logging
import functools
import warnings
import numpy
import sys
from coffea.analysis_tools import Weights

logger = logging.getLogger(__name__)


class DYStudiesProcessor(HggBaseProcessor):
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


class TagAndProbeProcessor(HggBaseProcessor):
    def __init__(
        self,
        metaconditions: Dict[str, Any],
        systematics: Dict[str, List[Any]] = None,
        corrections: Optional[Dict[str, List[str]]] = None,
        apply_trigger: bool = False,
        output_location: Optional[str] = None,
        taggers: Optional[List[Any]] = None,
        skipCQR: bool = False,
        skipJetVetoMap: bool = False,
        year: Optional[Dict[str, List[str]]] = None,
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
            trigger_group=".*SingleEle.*",
            analysis="tagAndProbe",
            skipCQR=skipCQR,
            skipJetVetoMap=False,
            year=year if year is not None else {},
            fiducialCuts=fiducialCuts,
            doDeco=doDeco,
            Smear_sigma_m=Smear_sigma_m,
            doFlow_corrections=doFlow_corrections,
            output_format=output_format
        )
        self.trigger_group = ".*SingleEle.*"
        self.analysis = "tagAndProbe"

        self.prefixes = {"tag": "tag", "probe": "probe"}

    def process(self, events: ak.Array) -> Dict[Any, Any]:

        dataset_name = events.metadata["dataset"]

        # data or mc?
        self.data_kind = "mc" if "GenPart" in ak.fields(events) else "data"

        # lumi mask
        if self.data_kind == "data":
            try:
                lumimask = select_lumis(self.year[dataset_name][0], events, logger)
                events = events[lumimask]
            except:
                logger.info(
                    f"[ lumimask ] Skip now! Unable to find year info of {dataset_name}"
                )
        # apply filters and triggers
        events = self.apply_filters_and_triggers(events)

        # we need ScEta for corrections and systematics, which is not present in NanoAODv11 but can be calculated using PV
        events.Photon = add_photon_SC_eta(events.Photon, events.PV)

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
            events = varying_function(events=events,year=self.year[dataset_name][0])

        for correction_name in correction_names:

            if correction_name in available_object_corrections.keys():
                logger.info(
                    f"\nApplying correction {correction_name} to dataset {dataset_name}\n"
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

        # Performing per photon corrections using normalizing flows
        if self.data_kind == "mc" and self.doFlow_corrections:

            # Applyting the Flow corrections to all photons before pre-selection
            counts = ak.num(original_photons)
            corrected_inputs,var_list = calculate_flow_corrections(original_photons, events, self.meta["flashggPhotons"]["flow_inputs"], self.meta["flashggPhotons"]["Isolation_transform_order"], year=self.year[dataset_name][0])

            # Store the raw nanoAOD value and update photon ID MVA value for preselection
            original_photons["mvaID_run3"] = ak.unflatten(self.add_photonid_mva_run3(original_photons, events), counts)
            original_photons["mvaID_nano"] = original_photons["mvaID"]

            # Store the raw values of the inputs and update the input values with the corrections since some variables used in the preselection
            for i in range(len(var_list)):
                original_photons["raw_" + str(var_list[i])] = original_photons[str(var_list[i])]
                original_photons[str(var_list[i])] = ak.unflatten(corrected_inputs[:,i] , counts)

            # Re-evaluate mvaID after corrections
            original_photons["mvaID"] = ak.unflatten(self.add_photonid_mva_run3(original_photons, events), counts)

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
                            systematic_dct["args"]["varying_function"], events=events, year=self.year[dataset_name][0]
                        )
                        # name=systematic_name, **systematic_dct["args"]
                    )
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
                photons_dct[f"{systematic}_{variation}"] = original_photons.systematics[
                    systematic
                ][variation]

        for variation, photons in photons_dct.items():
            logger.debug(f"Variation: {variation}")

            if variation == "nominal":
                do_variation = "nominal"

            if self.chained_quantile is not None:
                photons = self.chained_quantile.apply(photons, events)

            # recompute photonid_mva on the fly
            if self.photonid_mva_EB and self.photonid_mva_EE:
                photons = self.add_photonid_mva(photons, events)

            # photon preselection
            photons = photon_preselection(
                self, photons, events, apply_electron_veto=False, year=self.year[dataset_name][0]
            )

            if self.data_kind == "mc":
                # TODO: add weight systs and corrections! (if needed)
                # need to annotate the photons already here with a weight since later, each photon can be tag and probe and this changes the length of the array
                photons["weight"] = events["genWeight"]
                # keep only photons matched to gen e+ or e-
                photons = photons[photons.genPartFlav == 11]

                # make sure that the matched e+/e- comes from a Z
                gen_particles = events.GenPart
                gen_indices = photons.genPartIdx[photons.genPartIdx != -1]
                gen_particles = gen_particles[gen_indices]
                gen_particles = gen_particles[gen_particles.genPartIdxMother == 23]

            # other event related variables need to be added before the tag&probe combination
            # nPV just for validation of pileup reweighting
            photons["nPV"] = events.PV.npvs
            photons["fixedGridRhoAll"] = events.Rho.fixedGridRhoAll

            # TODO: HLT matching for data

            # double the number of diphoton candidates (each item in the pair can be both a tag and a probe)
            tnp = ak.combinations(photons, 2, fields=["tag", "probe"])
            pnt = ak.combinations(photons, 2, fields=["probe", "tag"])
            tnp_candidates = ak.concatenate([tnp, pnt], axis=1)

            # check that the e+/e- matched to tag and probe are not the same particle
            if self.data_kind == "mc":
                tnp_candidates = tnp_candidates[
                    tnp_candidates.tag.genPartIdx != tnp_candidates.probe.genPartIdx
                ]

            # tag selections
            tag_mask = (
                (tnp_candidates.tag.pt > 40)
                & (tnp_candidates.tag.electronIdx != -1)
                & (tnp_candidates.tag.pixelSeed)
                & (
                    tnp_candidates.tag.pfChargedIsoPFPV < 20
                )  # was: (tnp_candidates.tag.chargedHadronIso < 20)
                & (
                    tnp_candidates.tag.pfChargedIsoPFPV / tnp_candidates.tag.pt < 0.3
                )  # was: (tnp_candidates.tag.chargedHadronIso / tnp_candidates.tag.pt < 0.3)
            )

            # No selection on the probe to not bias it!

            # apply selections
            tnp_candidates = tnp_candidates[tag_mask]

            # Since the Weights object accepts only flat masks, the tag and probe mask is flattened
            flat_tag_and_probe_mask = ak.any(tag_mask, axis=1)

            """
            This n_event_tnp_cand array is created to keep track of how many tag and probe candidates we have at each event
            Since the pileup rw is calculated at a event level, we will have only one weight for event
            But since we are saving ak.flatten(tnp_candidates) , we need the n_event_tnp_cand to unroll the weights to each tnp candidate at the event
            """
            n_event_tnp_cand = [numpy.ones(n_tnp_candidates) for n_tnp_candidates in ak.num(tnp_candidates[flat_tag_and_probe_mask])]

            # candidates need to be flattened since we have each photon as a tag and probe, otherwise it can't be exported to numpy
            tnp_candidates = ak.flatten(tnp_candidates)

            # performing the weight corrections after the preselctions
            if self.data_kind == "mc":

                event_weights = Weights(size=len(events[flat_tag_and_probe_mask]))

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
                            events=events[flat_tag_and_probe_mask],
                            photons=events.Photon[flat_tag_and_probe_mask],
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
                            varying_function = available_weight_systematics[
                                systematic_name
                            ]
                            event_weights = varying_function(
                                events=events[flat_tag_and_probe_mask],
                                photons=events.Photon[flat_tag_and_probe_mask],
                                weights=event_weights,
                                dataset_name=dataset_name,
                                year=self.year[dataset_name][0],
                            )

            # Calculating sigma_m_overm_m
            if self.data_kind == "mc" and self.doFlow_corrections:
                tnp_candidates["sigma_m_over_m"] = 0.5 * numpy.sqrt(
                    (
                        tnp_candidates["tag"].raw_energyErr
                        / (
                            tnp_candidates["tag"].pt
                            * numpy.cosh(tnp_candidates["tag"].eta)
                        )
                    )
                    ** 2
                    + (
                        tnp_candidates["probe"].raw_energyErr
                        / (
                            tnp_candidates["probe"].pt
                            * numpy.cosh(tnp_candidates["probe"].eta)
                        )
                    )
                    ** 2
                )

                # We can now also calculate the corrected sigma_m_over_m
                tnp_candidates["sigma_m_over_m_corr"] = 0.5 * numpy.sqrt(
                    (
                        tnp_candidates["tag"].energyErr
                        / (
                            tnp_candidates["tag"].pt
                            * numpy.cosh(tnp_candidates["tag"].eta)
                        )
                    )
                    ** 2
                    + (
                        tnp_candidates["probe"].energyErr
                        / (
                            tnp_candidates["probe"].pt
                            * numpy.cosh(tnp_candidates["probe"].eta)
                        )
                    )
                    ** 2
                )

            else:
                tnp_candidates["sigma_m_over_m"] = 0.5 * numpy.sqrt(
                    (
                        tnp_candidates["tag"].energyErr
                        / (
                            tnp_candidates["tag"].pt
                            * numpy.cosh(tnp_candidates["tag"].eta)
                        )
                    )
                    ** 2
                    + (
                        tnp_candidates["probe"].energyErr
                        / (
                            tnp_candidates["probe"].pt
                            * numpy.cosh(tnp_candidates["probe"].eta)
                        )
                    )
                    ** 2
                )

            # Adding the smearing of the mass resolution also for the tag and probe workflow - for both data and Simulation!
            # Just a reminder, the pt/energy of teh data is not smearing, but the smearing term is added to the data sigma_m_over_m
            if self.Smear_sigma_m:

                # Adding the smeared energyErr error to the ntuples!
                if self.doFlow_corrections and self.data_kind == "mc":
                    tnp_candidates["tag","energyErr_Smeared"] = numpy.sqrt((tnp_candidates["tag"].raw_energyErr)**2 + (tnp_candidates["tag"].rho_smear * ((tnp_candidates["tag"].pt * numpy.cosh(tnp_candidates["tag"].eta)))) ** 2)
                    tnp_candidates["probe","energyErr_Smeared"] = numpy.sqrt((tnp_candidates["probe"].raw_energyErr) ** 2 + (tnp_candidates["probe"].rho_smear * ((tnp_candidates["probe"].pt * numpy.cosh(tnp_candidates["probe"].eta)))) ** 2)

                    tnp_candidates["sigma_m_over_m_Smeared"] = 0.5 * numpy.sqrt(
                        (
                            tnp_candidates["tag"].energyErr_Smeared / (tnp_candidates["tag"].pt * numpy.cosh(tnp_candidates["tag"].eta))
                        )
                        ** 2
                        + (
                            tnp_candidates["probe"].energyErr_Smeared / (tnp_candidates["probe"].pt * numpy.cosh(tnp_candidates["probe"].eta))
                        )
                        ** 2
                    )

                    # Adding the smeared energyErr error to the ntuples!
                    tnp_candidates["tag","corr_energyErr_Smeared"] = numpy.sqrt((tnp_candidates["tag"].energyErr) ** 2 + (tnp_candidates["tag"].rho_smear * ((tnp_candidates["tag"].pt * numpy.cosh(tnp_candidates["tag"].eta)))) ** 2)
                    tnp_candidates["probe","corr_energyErr_Smeared"] = numpy.sqrt((tnp_candidates["probe"].energyErr) ** 2 + (tnp_candidates["probe"].rho_smear * ((tnp_candidates["probe"].pt * numpy.cosh(tnp_candidates["probe"].eta)))) ** 2)

                    tnp_candidates["sigma_m_over_m_Smeared_corrected"] = 0.5 * numpy.sqrt(
                        (
                            tnp_candidates["tag"].corr_energyErr_Smeared / (tnp_candidates["tag"].pt * numpy.cosh(tnp_candidates["tag"].eta))
                        )
                        ** 2
                        + (
                            tnp_candidates["probe"].corr_energyErr_Smeared / (tnp_candidates["probe"].pt * numpy.cosh(tnp_candidates["probe"].eta))
                        )
                        ** 2
                    )

                else:
                    tnp_candidates["tag","energyErr_Smeared"] = numpy.sqrt((tnp_candidates["tag"].energyErr)**2 + (tnp_candidates["tag"].rho_smear * ((tnp_candidates["tag"].pt * numpy.cosh(tnp_candidates["tag"].eta)))) ** 2)
                    tnp_candidates["probe","energyErr_Smeared"] = numpy.sqrt((tnp_candidates["probe"].energyErr) ** 2 + (tnp_candidates["probe"].rho_smear * ((tnp_candidates["probe"].pt * numpy.cosh(tnp_candidates["probe"].eta)))) ** 2)

                    tnp_candidates["sigma_m_over_m_Smeared"] = 0.5 * numpy.sqrt(
                        (
                            tnp_candidates["tag"].energyErr_Smeared / (tnp_candidates["tag"].pt * numpy.cosh(tnp_candidates["tag"].eta))
                        )
                        ** 2
                        + (
                            tnp_candidates["probe"].energyErr_Smeared / (tnp_candidates["probe"].pt * numpy.cosh(tnp_candidates["probe"].eta))
                        )
                        ** 2
                    )

            # Adding the tagandprobe pair mass. Based on the expression provided here for a massless pair of particles -> (https://en.wikipedia.org/wiki/Invariant_mass)
            tnp_candidates["mass"] = numpy.sqrt(2 * tnp_candidates["tag"].pt * tnp_candidates["probe"].pt * (numpy.cosh(tnp_candidates["tag"].eta - tnp_candidates["probe"].eta) - numpy.cos(tnp_candidates["tag"].phi - tnp_candidates["probe"].phi)))

            if self.output_location is not None:
                df = diphoton_list_to_pandas(self, tnp_candidates)

                # since we annotated the photons with event variables, these exist now for tag and probe. This concerns weights as well as nPV and fixedGridRhoAll Remove:
                if self.data_kind == "mc":

                    # Store variations with respect to central weight
                    if do_variation == "nominal":
                        if len(event_weights.variations):
                            logger.info(
                                "Adding systematic weight variations to nominal output file."
                            )
                        for modifier in event_weights.variations:
                            df["weight_" + modifier] = numpy.hstack(event_weights.weight(
                                modifier=modifier
                            ) * n_event_tnp_cand)

                    # storing the central weights
                    df["weight_central"] = numpy.hstack(event_weights.weight() * n_event_tnp_cand)
                    # generated weights * other weights (pile up, SF, etc ...)
                    df["weight"] = df["tag_weight"] * numpy.hstack(event_weights.weight() * n_event_tnp_cand)
                    df["weight_no_pu"] = df["tag_weight"]

                    # dropping the nominal and varitation weights
                    df.drop(["tag_weight", "probe_weight"], axis=1, inplace=True)

                df["nPV"] = df["tag_nPV"]
                df.drop(["tag_nPV", "probe_nPV"], axis=1, inplace=True)
                df["fixedGridRhoAll"] = df["tag_fixedGridRhoAll"]
                df.drop(
                    ["tag_fixedGridRhoAll", "probe_fixedGridRhoAll"],
                    axis=1,
                    inplace=True,
                )

                fname = (
                    events.behavior["__events_factory__"]._partition_key.replace(
                        "/", "_"
                    )
                    + ".%s" % self.output_format
                )
                subdirs = []
                if "dataset" in events.metadata:
                    subdirs.append(events.metadata["dataset"])
                subdirs.append(variation)
                dump_pandas(self, df, fname, self.output_location, subdirs)

        return {}

    def process_extra(self, events: ak.Array) -> ak.Array:
        return events, {}

    def postprocess(self, accumulant: Dict[Any, Any]) -> Any:
        pass
