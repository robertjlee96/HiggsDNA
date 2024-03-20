from higgs_dna.workflows.base import HggBaseProcessor
from higgs_dna.systematics import object_corrections as available_object_corrections
from higgs_dna.systematics import weight_corrections as available_weight_corrections
from higgs_dna.utils.dumping_utils import diphoton_ak_array, dump_ak_array, diphoton_list_to_pandas, dump_pandas

from typing import Any, Dict, List, Optional
import awkward
import logging
import warnings
import numpy
import sys
from coffea.analysis_tools import Weights

logger = logging.getLogger(__name__)


class ParticleLevelProcessor(HggBaseProcessor):
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
        fiducialCuts: str = "none",
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

    def process_extra(self, events: awkward.Array) -> awkward.Array:
        return events, {}

    def process(self, events: awkward.Array) -> Dict[Any, Any]:
        dataset_name = events.metadata["dataset"]

        # metadata array to append to higgsdna output
        metadata = {}

        # data or monte carlo?
        self.data_kind = "mc" if hasattr(events, "GenPart") else "data"

        if self.data_kind == "data":
            logger.info("The 'particleLevel' processor can only be run on MC. Aborting now...")
            sys.exit(0)

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

        # Add sum of gen weights before selection for normalisation in postprocessing
        metadata["sum_genw_presel"] = str(awkward.sum(events.genWeight))

        # read which systematics and corrections to process
        try:
            correction_names = self.corrections[dataset_name]
        except KeyError:
            correction_names = []

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

        # Even though it is a particle-level processor, we will at least save basic diphoton quantities on reco level for comparison
        photons = events.Photon
        photons = photons[awkward.argsort(photons.pt, ascending=False)]
        photons["charge"] = awkward.zeros_like(
            photons.pt
        )
        diphotons = awkward.combinations(
            photons, 2, fields=["pho_lead", "pho_sublead"]
        )

        # now turn the diphotons into candidates with four momenta and such
        diphoton_4mom = diphotons["pho_lead"] + diphotons["pho_sublead"]
        diphotons["pt"] = diphoton_4mom.pt
        diphotons["eta"] = diphoton_4mom.eta
        diphotons["phi"] = diphoton_4mom.phi
        diphotons["mass"] = diphoton_4mom.mass
        diphotons["charge"] = diphoton_4mom.charge

        diphoton_pz = diphoton_4mom.z
        diphoton_e = diphoton_4mom.energy
        diphotons["rapidity"] = 0.5 * numpy.log((diphoton_e + diphoton_pz) / (diphoton_e - diphoton_pz))

        diphotons = awkward.with_name(diphotons, "PtEtaPhiMCandidate")

        # sort diphotons by pT
        diphotons = diphotons[
            awkward.argsort(diphotons.pt, ascending=False)
        ]

        events["GenIsolatedPhoton"] = awkward.pad_none(events["GenIsolatedPhoton"], 2)

        diphotons["leadingGenIsolatedPhoton_pt"] = events.GenIsolatedPhoton[:,0].pt
        diphotons["leadingGenIsolatedPhoton_eta"] = events.GenIsolatedPhoton[:,0].eta
        diphotons["leadingGenIsolatedPhoton_phi"] = events.GenIsolatedPhoton[:,0].phi

        diphotons["subleadingGenIsolatedPhoton_pt"] = events.GenIsolatedPhoton[:,1].pt
        diphotons["subleadingGenIsolatedPhoton_eta"] = events.GenIsolatedPhoton[:,1].eta
        diphotons["subleadingGenIsolatedPhoton_phi"] = events.GenIsolatedPhoton[:,1].phi

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

        # workflow specific processing
        events, process_extra = self.process_extra(events)
        histos_etc.update(process_extra)

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
        events["diphotons"] = diphotons
        # annotate diphotons with event information
        diphotons["event"] = events.event
        diphotons["lumi"] = events.luminosityBlock
        diphotons["run"] = events.run
        # annotate diphotons with dZ information (difference between z position of GenVtx and PV) as required by flashggfinalfits
        diphotons["genWeight"] = events.genWeight
        diphotons["dZ"] = events.GenVtx.z - events.PV.z
        # Necessary for differential xsec measurements in final fits ("truth" variables)
        diphotons["HTXS_Higgs_pt"] = events.HTXS.Higgs_pt
        diphotons["HTXS_Higgs_y"] = events.HTXS.Higgs_y
        diphotons["HTXS_njets30"] = events.HTXS.njets30  # Need to clarify if this variable is suitable, does it fulfill abs(eta_j) < 2.5? Probably not
        # Preparation for HTXS measurements later, start with stage 0 to disentangle VH into WH and ZH for final fits
        diphotons["HTXS_stage_0"] = events.HTXS.stage_0

        # return if there is no surviving events
        if len(diphotons) == 0:
            logger.debug("No surviving events in this run, return now!")
            return histos_etc

        # Retain all events
        selection_mask = diphotons.pho_lead.pt > -10
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
                    photons=events["diphotons"][
                        selection_mask
                    ],
                    weights=event_weights,
                    dataset_name=dataset_name,
                    year=self.year[dataset_name][0],
                )
        diphotons["weight_central"] = event_weights.weight()

        # Multiply weight by genWeight for normalisation in post-processing chain
        event_weights._weight = (
            events["genWeight"][selection_mask]
            * diphotons["weight_central"]
        )
        diphotons["weight"] = event_weights.weight()

        if self.output_location is not None:
            if self.output_format == "root":
                df = diphoton_list_to_pandas(self, diphotons)
            else:
                akarr = diphoton_ak_array(self, diphotons)
            fname = (
                events.behavior[
                    "__events_factory__"
                ]._partition_key.replace("/", "_")
                + ".%s" % self.output_format
            )
            subdirs = []
            if "dataset" in events.metadata:
                subdirs.append(events.metadata["dataset"])
            if self.output_format == "root":
                dump_pandas(self, df, fname, self.output_location, subdirs)
            else:
                dump_ak_array(
                    self, akarr, fname, self.output_location, metadata, subdirs,
                )

        return histos_etc

    def postprocess(self, accumulant: Dict[Any, Any]) -> Any:
        pass
