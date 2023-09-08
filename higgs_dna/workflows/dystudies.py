from higgs_dna.workflows.base import HggBaseProcessor
from higgs_dna.systematics import object_systematics as available_object_systematics
from higgs_dna.systematics import object_corrections as available_object_corrections
from higgs_dna.systematics import weight_systematics as available_weight_systematics
from higgs_dna.systematics import weight_corrections as available_weight_corrections
from higgs_dna.selections.photon_selections import photon_preselection
from higgs_dna.utils.dumping_utils import diphoton_list_to_pandas, dump_pandas
from higgs_dna.tools.SC_eta import add_photon_SC_eta
from typing import Any, Dict, List, Optional
import awkward as ak
import logging
import functools
import warnings
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
            output_format=output_format,
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
        output_format: str = "parquet",
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
            output_format=output_format,
        )
        self.trigger_group = ".*SingleEle.*"
        self.analysis = "tagAndProbe"

        self.prefixes = {"tag": "tag", "probe": "probe"}

    def process(self, events: ak.Array) -> Dict[Any, Any]:

        dataset_name = events.metadata["dataset"]

        # data or mc?
        self.data_kind = "mc" if "GenPart" in ak.fields(events) else "data"

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

        # object corrections:
        for correction_name in correction_names:
            if correction_name in available_object_corrections.keys():
                logger.info(
                    f"\nApplying correction {correction_name} to dataset {dataset_name}\n"
                )
                varying_function = available_object_corrections[correction_name]
                events = varying_function(events=events)
            elif correction_name in available_weight_corrections:
                # event weight corrections will be applied after photon preselection / application of further taggers
                continue
            else:
                # may want to throw an error instead, needs to be discussed
                warnings.warn(f"Could not process correction {correction_name}.")
                continue

        original_photons = events.Photon
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

        if self.data_kind == "mc":
            event_weights = Weights(size=len(events))
            # _weight will correspond to "nominal" weight, what else has to be included here? (lumi? xSec? MC sum of weights?)
            event_weights._weight = events["genWeight"]

        for variation, photons in photons_dct.items():
            logger.debug(f"Variation: {variation}")

            if self.chained_quantile is not None:
                photons = self.chained_quantile.apply(photons, events)

            # recompute photonid_mva on the fly
            if self.photonid_mva_EB and self.photonid_mva_EE:
                photons = self.add_photonid_mva(photons, events)

            # photon preselection
            photons = photon_preselection(
                self, photons, events, apply_electron_veto=False
            )

            if self.data_kind == "mc":
                # TODO: add weight systs and corrections! (if needed)
                # need to annotate the photons already here with a weight since later, each photon can be tag and probe and this changes the length of the array
                photons["weight"] = event_weights.weight()
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

            # probe selections
            probe_mask = (
                tnp_candidates.probe.pfChargedIsoPFPV < 20
            ) & (  # was: (tnp_candidates.probe.chargedHadronIso < 20)
                tnp_candidates.probe.pfChargedIsoPFPV / tnp_candidates.probe.pt
                < 0.3  # was: tnp_candidates.probe.chargedHadronIso / tnp_candidates.probe.pt < 0.3
            )

            # apply selections
            tnp_candidates = tnp_candidates[tag_mask & probe_mask]
            # candidates need to be flattened since we have each photon as a tag and probe, otherwise it can't be exported to numpy
            tnp_candidates = ak.flatten(tnp_candidates)

            if self.output_location is not None:
                df = diphoton_list_to_pandas(self, tnp_candidates)

                # since we annotated the photons with event variables, these exist now for tag and probe. This concerns weights as well as nPV and fixedGridRhoAll Remove:
                if self.data_kind == "mc":
                    df["weight"] = df["tag_weight"]
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
