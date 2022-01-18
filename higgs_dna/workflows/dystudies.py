from higgs_dna.workflows.base import HggBaseProcessor

from typing import Any, Dict, List, Optional
import awkward as ak
import logging

logger = logging.getLogger(__name__)


class DYStudiesProcessor(HggBaseProcessor):
    def __init__(
        self,
        metaconditions: Dict[str, Any],
        do_systematics: bool = False,
        apply_trigger: bool = False,
        output_location: Optional[str] = None,
        taggers: Optional[List[Any]] = None,
        skipCQR: bool = False,
    ) -> None:
        super().__init__(
            metaconditions,
            do_systematics=do_systematics,
            apply_trigger=apply_trigger,
            output_location=output_location,
            taggers=taggers,
            trigger_group=".*DoubleEG.*",
            analysis="mainAnalysis",
            skipCQR=skipCQR,
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
        do_systematics: bool = False,
        apply_trigger: bool = False,
        output_location: Optional[str] = None,
        taggers: Optional[List[Any]] = None,
        skipCQR: bool = False,
    ) -> None:
        super().__init__(
            metaconditions,
            do_systematics=do_systematics,
            apply_trigger=apply_trigger,
            output_location=output_location,
            taggers=taggers,
            trigger_group=".*SingleEle.*",
            analysis="tagAndProbe",
            skipCQR=skipCQR,
        )
        self.trigger_group = ".*SingleEle.*"
        self.analysis = "tagAndProbe"

        self.prefixes = {"tag": "tag", "probe": "probe"}

    def process(self, events: ak.Array) -> Dict[Any, Any]:
        # data or mc?
        self.data_kind = "mc" if "GenPart" in ak.fields(events) else "data"

        # apply filters and triggers
        events = self.apply_filters_and_triggers(events)

        # modifications to photons
        photons = events.Photon

        if self.chained_quantile is not None:
            photons = self.chained_quantile.apply(events)

        # recompute photonid_mva on the fly
        if self.photonid_mva_EB and self.photonid_mva_EE:
            photons = self.add_photonid_mva(photons, events)

        # photon preselection
        photons = self.photon_preselection(photons, events)

        if self.data_kind == "mc":
            # keep only photons matched to gen e+ or e-
            photons = photons[photons.genPartFlav == 11]

            # make sure that the matched e+/e- comes from a Z
            gen_particles = events.GenPart
            gen_indices = photons.genPartIdx[photons.genPartIdx != -1]
            gen_particles = gen_particles[gen_indices]
            gen_particles = gen_particles[gen_particles.genPartIdxMother == 23]

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
            & (tnp_candidates.tag.chargedHadronIso < 20)
            & (tnp_candidates.tag.chargedHadronIso / tnp_candidates.tag.pt < 0.3)
        )

        # probe selections
        probe_mask = (tnp_candidates.probe.chargedHadronIso < 20) & (
            tnp_candidates.probe.chargedHadronIso / tnp_candidates.probe.pt < 0.3
        )

        # apply selections
        tnp_candidates = tnp_candidates[tag_mask & probe_mask]

        # candidates need to be flattened since we have each photon as a tag and probe, otherwise it can't be exported to numpy
        tnp_candidates = ak.flatten(tnp_candidates)

        if self.output_location is not None:
            df = self.diphoton_list_to_pandas(tnp_candidates)
            fname = (
                events.behavior["__events_factory__"]._partition_key.replace("/", "_")
                + ".parquet"
            )
            subdirs = []
            if "dataset" in events.metadata:
                subdirs.append(events.metadata["dataset"])
            self.dump_pandas(df, fname, self.output_location, subdirs)

        return {}

    def process_extra(self, events: ak.Array) -> ak.Array:
        return events, {}

    def postprocess(self, accumulant: Dict[Any, Any]) -> Any:
        pass
