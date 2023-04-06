from higgs_dna.tools.chained_quantile import ChainedQuantileRegression
from higgs_dna.tools.diphoton_mva import calculate_diphoton_mva
from higgs_dna.tools.xgb_loader import load_bdt
from higgs_dna.tools.photonid_mva import calculate_photonid_mva, load_photonid_mva
from higgs_dna.metaconditions import photon_id_mva_weights
from higgs_dna.metaconditions import diphoton as diphoton_mva_dir
from higgs_dna.systematics import object_systematics as available_object_systematics
from higgs_dna.systematics import object_corrections as available_object_corrections
from higgs_dna.systematics import weight_systematics as available_weight_systematics
from higgs_dna.systematics import weight_corrections as available_weight_corrections

import functools
import operator
import os
import pathlib
import shutil
import warnings
from typing import Any, Dict, List, Optional
import awkward
import numpy
import pandas
import vector
from coffea import processor
from coffea.analysis_tools import Weights

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
    ) -> None:
        self.meta = metaconditions
        self.systematics = systematics if systematics is not None else {}
        self.corrections = corrections if corrections is not None else {}
        self.apply_trigger = apply_trigger
        self.output_location = output_location
        self.trigger_group = trigger_group
        self.analysis = analysis
        self.skipCQR = skipCQR

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
        self.max_trkSumPtHollowConeDR03_EB_low_r9 = 6.0  # for v11, we cut on Photon_pfChargedIsoPFPV
        self.max_trkSumPtHollowConeDR03_EE_low_r9 = 6.0  # Leaving the names of the preselection cut variables the same to change as little as possible
        self.max_sieie_EB_low_r9 = 0.015
        self.max_sieie_EE_low_r9 = 0.035
        self.max_pho_iso_EB_low_r9 = 4.0
        self.max_pho_iso_EE_low_r9 = 4.0

        self.eta_rho_corr = 1.5
        self.low_eta_rho_corr = 0.16544
        self.high_eta_rho_corr = 0.13212

        logger.debug(f"Setting up processor with metaconditions: {self.meta}")

        self.taggers = []
        if taggers is not None:
            self.taggers = taggers
            self.taggers.sort(key=lambda x: x.priority)

        self.prefixes = {"pho_lead": "lead", "pho_sublead": "sublead"}

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

    def photon_preselection(
        self, photons: awkward.Array, events: awkward.Array
    ) -> awkward.Array:
        """
        Apply preselection cuts to photons.
        Note that these selections are applied on each photon, it is not based on the diphoton pair.
        """
        # hlt-mimicking cuts
        rho = events.Rho.fixedGridRhoAll * awkward.ones_like(photons.pt)
        photon_abs_eta = numpy.abs(photons.eta)
        isEB_high_r9 = (photon_abs_eta < self.gap_barrel_eta) & (
            photons.r9 > self.min_full5x5_r9_EB_high_r9
        )
        isEE_high_r9 = (photon_abs_eta > self.gap_endcap_eta) & (
            photons.r9 > self.min_full5x5_r9_EE_high_r9
        )
        isEB_low_r9 = (
            (photon_abs_eta < self.gap_barrel_eta)
            & (photons.r9 > self.min_full5x5_r9_EB_low_r9)
            & (photons.r9 < self.min_full5x5_r9_EB_high_r9)
            & (
                photons.pfChargedIsoPFPV  # for v11
                < self.max_trkSumPtHollowConeDR03_EB_low_r9
            )
            & (photons.sieie < self.max_sieie_EB_low_r9)
            & (
                (
                    (photon_abs_eta < self.eta_rho_corr)
                    & (
                        photons.pfPhoIso03 - rho * self.low_eta_rho_corr
                        < self.max_pho_iso_EB_low_r9
                    )
                )
                | (
                    (photon_abs_eta > self.eta_rho_corr)
                    & (
                        photons.pfPhoIso03 - rho * self.high_eta_rho_corr
                        < self.max_pho_iso_EB_low_r9
                    )
                )
            )
        )
        isEE_low_r9 = (
            (photon_abs_eta < self.gap_barrel_eta)
            & (photons.r9 > self.min_full5x5_r9_EE_low_r9)
            & (photons.r9 < self.min_full5x5_r9_EE_high_r9)
            & (
                photons.pfChargedIsoPFPV  # for v11
                < self.max_trkSumPtHollowConeDR03_EE_low_r9
            )
            & (photons.sieie < self.max_sieie_EE_low_r9)
            & (
                (
                    (photon_abs_eta < self.eta_rho_corr)
                    & (
                        photons.pfPhoIso03 - rho * self.low_eta_rho_corr
                        < self.max_pho_iso_EE_low_r9
                    )
                )
                | (
                    (photon_abs_eta > self.eta_rho_corr)
                    & (
                        photons.pfPhoIso03 - rho * self.high_eta_rho_corr
                        < self.max_pho_iso_EE_low_r9
                    )
                )
            )
        )

        return photons[
            (photons.pt > self.min_pt_photon)
            & (photon_abs_eta < self.max_sc_eta)
            & (
                (photon_abs_eta < self.gap_barrel_eta)
                | (photon_abs_eta > self.gap_endcap_eta)
            )
            & (photons.mvaID > self.min_mvaid)
            & (photons.hoe < self.max_hovere)
            & (
                (photons.r9 > self.min_full5x5_r9)
                | (photons.pfRelIso03_chg_quadratic < self.max_chad_iso)  # changed from pfRelIso03_chg since this variable is not in v11 nanoAOD...?
                | (photons.pfRelIso03_chg_quadratic / photons.pt < self.max_chad_rel_iso)
            )
            & (isEB_high_r9 | isEB_low_r9 | isEE_high_r9 | isEE_low_r9)
        ]

    def diphoton_list_to_pandas(self, diphotons: awkward.Array) -> pandas.DataFrame:
        """
        Convert diphoton array to pandas dataframe.
        By default the observables related to each item of the diphoton pair are
        stored preceded by its prefix (e.g. 'lead', 'sublead').
        The observables related to the diphoton pair are stored with no prefix.
        To change the behavior, you can redefine the `diphoton_list_to_pandas` method in the
        derived class.
        """
        output = pandas.DataFrame()
        for field in awkward.fields(diphotons):
            prefix = self.prefixes.get(field, "")
            if len(prefix) > 0:
                for subfield in awkward.fields(diphotons[field]):
                    if subfield != "__systematics__":
                        output[f"{prefix}_{subfield}"] = awkward.to_numpy(
                            diphotons[field][subfield]
                        )
            else:
                output[field] = awkward.to_numpy(diphotons[field])
        return output

    def dump_pandas(
        self,
        pddf: pandas.DataFrame,
        fname: str,
        location: str,
        subdirs: Optional[List[str]] = None,
    ) -> None:
        """
        Dump a pandas dataframe to disk at location/'/'.join(subdirs)/fname.
        """
        subdirs = subdirs or []
        xrd_prefix = "root://"
        pfx_len = len(xrd_prefix)
        xrootd = False
        if xrd_prefix in location:
            try:
                import XRootD  # type: ignore
                import XRootD.client  # type: ignore

                xrootd = True
            except ImportError as err:
                raise ImportError(
                    "Install XRootD python bindings with: conda install -c conda-forge xroot"
                ) from err
        local_file = (
            os.path.abspath(os.path.join(".", fname))
            if xrootd
            else os.path.join(".", fname)
        )
        merged_subdirs = "/".join(subdirs) if xrootd else os.path.sep.join(subdirs)
        destination = (
            location + merged_subdirs + f"/{fname}"
            if xrootd
            else os.path.join(location, os.path.join(merged_subdirs, fname))
        )
        pddf.to_parquet(local_file)
        if xrootd:
            copyproc = XRootD.client.CopyProcess()
            copyproc.add_job(local_file, destination)
            copyproc.prepare()
            copyproc.run()
            client = XRootD.client.FileSystem(
                location[: location[pfx_len:].find("/") + pfx_len]
            )
            status = client.locate(
                destination[destination[pfx_len:].find("/") + pfx_len + 1 :],
                XRootD.client.flags.OpenFlags.READ,
            )
            assert status[0].ok
            del client
            del copyproc
        else:
            dirname = os.path.dirname(destination)
            if not os.path.exists(dirname):
                pathlib.Path(dirname).mkdir(parents=True, exist_ok=True)
            shutil.copy(local_file, destination)
            assert os.path.isfile(destination)
        pathlib.Path(local_file).unlink()

    def diphoton_ak_array(self, diphotons: awkward.Array) -> awkward.Array:
        """
        Adjust the prefix.
        By default the observables related to each item of the diphoton pair are
        stored preceded by its prefix (e.g. 'lead', 'sublead').
        The observables related to the diphoton pair are stored with no prefix.
        """
        output = {}
        for field in awkward.fields(diphotons):
            prefix = self.prefixes.get(field, "")
            if len(prefix) > 0:
                for subfield in awkward.fields(diphotons[field]):
                    if subfield != "__systematics__":
                        output[f"{prefix}_{subfield}"] = diphotons[field][subfield]
            else:
                output[field] = diphotons[field]
        return awkward.Array(output)

    def dump_ak_array(
        self,
        akarr: awkward.Array,
        fname: str,
        location: str,
        subdirs: Optional[List[str]] = None,
    ) -> None:
        """
        Dump a pandas dataframe to disk at location/'/'.join(subdirs)/fname.
        """
        subdirs = subdirs or []
        xrd_prefix = "root://"
        pfx_len = len(xrd_prefix)
        xrootd = False
        if xrd_prefix in location:
            try:
                import XRootD  # type: ignore
                import XRootD.client  # type: ignore

                xrootd = True
            except ImportError as err:
                raise ImportError(
                    "Install XRootD python bindings with: conda install -c conda-forge xroot"
                ) from err
        local_file = (
            os.path.abspath(os.path.join(".", fname))
            if xrootd
            else os.path.join(".", fname)
        )
        merged_subdirs = "/".join(subdirs) if xrootd else os.path.sep.join(subdirs)
        destination = (
            location + merged_subdirs + f"/{fname}"
            if xrootd
            else os.path.join(location, os.path.join(merged_subdirs, fname))
        )
        awkward.to_parquet(akarr, local_file)
        if xrootd:
            copyproc = XRootD.client.CopyProcess()
            copyproc.add_job(local_file, destination)
            copyproc.prepare()
            copyproc.run()
            client = XRootD.client.FileSystem(
                location[: location[pfx_len:].find("/") + pfx_len]
            )
            status = client.locate(
                destination[destination[pfx_len:].find("/") + pfx_len + 1 :],
                XRootD.client.flags.OpenFlags.READ,
            )
            assert status[0].ok
            del client
            del copyproc
        else:
            dirname = os.path.dirname(destination)
            pathlib.Path(dirname).mkdir(parents=True, exist_ok=True)
            shutil.copy(local_file, destination)
            assert os.path.isfile(destination)
        pathlib.Path(local_file).unlink()

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
        # data or monte carlo?
        self.data_kind = "mc" if hasattr(events, "GenPart") else "data"

        # apply filters and triggers
        events = self.apply_filters_and_triggers(events)

        # here we start recording possible coffea accumulators
        # most likely histograms, could be counters, arrays, ...
        histos_etc = {}

        # read which systematics and corrections to process
        dataset_name = events.metadata["dataset"]
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
                        name=systematic_name, kind=systematic_dct["args"]["kind"], what=systematic_dct["args"]["what"], varying_function=functools.partial(systematic_dct["args"]["varying_function"], events=events)
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
                photons_dct[f"{systematic}_{variation}"] = original_photons.systematics[
                    systematic
                ][variation]

        for variation, photons in photons_dct.items():
            logger.debug(f"Variation: {variation}")
            if self.chained_quantile is not None:
                photons = self.chained_quantile.apply(photons, events)

            # recompute photonid_mva on the fly
            if self.photonid_mva_EB and self.photonid_mva_EE:
                photons = self.add_photonid_mva(photons, events)

            # photon preselection
            photons = self.photon_preselection(photons, events)
            # sort photons in each event descending in pt
            # make descending-pt combinations of photons
            photons = photons[awkward.argsort(photons.pt, ascending=False)]
            photons["charge"] = awkward.zeros_like(photons.pt)  # added this because charge is not a property of photons in nanoAOD v11. We just assume every photon has charge zero...
            diphotons = awkward.combinations(
                photons, 2, fields=["pho_lead", "pho_sublead"]
            )
            # the remaining cut is to select the leading photons
            # the previous sort assures the order
            diphotons = diphotons[diphotons["pho_lead"].pt > self.min_pt_lead_photon]

            # now turn the diphotons into candidates with four momenta and such
            diphoton_4mom = diphotons["pho_lead"] + diphotons["pho_sublead"]
            diphotons["pt"] = diphoton_4mom.pt
            diphotons["eta"] = diphoton_4mom.eta
            diphotons["phi"] = diphoton_4mom.phi
            diphotons["mass"] = diphoton_4mom.mass
            diphotons["charge"] = diphoton_4mom.charge
            diphotons = awkward.with_name(diphotons, "PtEtaPhiMCandidate")

            # sort diphotons by pT
            diphotons = diphotons[awkward.argsort(diphotons.pt, ascending=False)]

            # baseline modifications to diphotons
            if self.diphoton_mva is not None:
                diphotons = self.add_diphoton_mva(diphotons, events)

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
                            diphotons["_".join([tagger.name, str(tagger.priority)])]
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
            events[f"diphotons_{variation}"] = diphotons
            # event PDF/Scale/PS weights
            if hasattr(events, "LHEScaleWeight"):
                diphotons["nLHEScaleWeight"] = awkward.num(events.LHEScaleWeight,axis=1)
                diphotons["LHEScaleWeight"] = events.LHEScaleWeight
            if hasattr(events, "LHEPdfWeight"):
                diphotons["nLHEPdfWeight"] = awkward.num(events.LHEPdfWeight,axis=1)
                diphotons["LHEPdfWeight"] = events.LHEPdfWeight
            if hasattr(events, "PSWeight"):
                diphotons["nPSWeight"] = awkward.num(events.PSWeight,axis=1)
                diphotons["PSWeight"] = events.PSWeight

            # annotate diphotons with event information
            diphotons["event"] = events.event
            diphotons["lumi"] = events.luminosityBlock
            diphotons["run"] = events.run
            # annotate diphotons with dZ information (difference between z position of GenVtx and PV) as required by flashggfinalfits
            if self.data_kind == "mc":
                diphotons["dZ"] = events.GenVtx.z - events.PV.z
            # Fill zeros for data because there is no GenVtx for data, obviously
            else:
                diphotons["dZ"] = awkward.zeros_like(events.PV.z)

            # drop events without a preselected diphoton candidate
            # drop events without a tag, if there are tags
            if len(self.taggers):
                selection_mask = ~(
                    awkward.is_none(diphotons) | awkward.is_none(diphotons.best_tag)
                )
                diphotons = diphotons[selection_mask]
            else:
                selection_mask = ~awkward.is_none(diphotons)
                diphotons = diphotons[selection_mask]

            if self.data_kind == "mc":
                # initiate Weight container here, after selection, since event selection cannot easily be applied to weight container afterwards
                event_weights = Weights(size=len(events[selection_mask]))
                # _weight will correspond to "nominal" weight, what else has to be included here? (lumi? xSec? MC sum of weights?)
                event_weights._weight = events["genWeight"][selection_mask]

                # corrections to event weights:
                for correction_name in correction_names:
                    if correction_name in available_weight_corrections:
                        logger.info(
                            f"Adding correction {correction_name} to weight collection of dataset {dataset_name}"
                        )
                        varying_function = available_weight_corrections[correction_name]
                        event_weights = varying_function(
                            events=events[selection_mask],
                            photons=events[f"diphotons_{variation}"][selection_mask],
                            weights=event_weights,
                        )

                # systematic variations of event weights go to nominal output dataframe:
                if variation == "nominal":
                    for systematic_name in systematic_names:
                        if systematic_name in available_weight_systematics:
                            logger.info(
                                f"Adding systematic {systematic_name} to weight collection of dataset {dataset_name}"
                            )
                            varying_function = available_weight_systematics[
                                systematic_name
                            ]
                            event_weights = varying_function(
                                events=events[selection_mask],
                                photons=events[f"diphotons_{variation}"][
                                    selection_mask
                                ],
                                weights=event_weights,
                            )

                diphotons["weight"] = event_weights.weight()
                if variation == "nominal":
                    if len(event_weights.variations):
                        logger.info(
                            "Adding systematic weight variations to nominal output file."
                        )
                    for modifier in event_weights.variations:
                        diphotons["weight_" + modifier] = event_weights.weight(
                            modifier=modifier
                        )

            if self.output_location is not None:
                # df = self.diphoton_list_to_pandas(diphotons)
                akarr = self.diphoton_ak_array(diphotons)
                fname = (
                    events.behavior["__events_factory__"]._partition_key.replace(
                        "/", "_"
                    )
                    + ".parquet"
                )
                subdirs = []
                if "dataset" in events.metadata:
                    subdirs.append(events.metadata["dataset"])
                subdirs.append(variation)
                # self.dump_pandas(df, fname, self.output_location, subdirs)
                self.dump_ak_array(akarr, fname, self.output_location, subdirs)

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