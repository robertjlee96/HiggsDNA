from higgs_dna.tools.chained_quantile import ChainedQuantileRegression
from higgs_dna.tools.diphoton_mva import calculate_diphoton_mva
from higgs_dna.tools.xgb_loader import load_bdt
from higgs_dna.tools.photonid_mva import calculate_photonid_mva, load_photonid_mva
from higgs_dna.metaconditions import photon_id_mva_weights
from higgs_dna.metaconditions import diphoton as diphoton_mva_dir

from dataclasses import dataclass
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

import logging

logger = logging.getLogger(__name__)

vector.register_awkward()


@dataclass
class HggBaseProcessor(processor.ProcessorABC):  # type: ignore
    meta: Dict[str, Any]
    do_systematics: bool
    apply_trigger: bool
    trigger_group: str
    analysis: str
    output_location: Optional[str] = None
    taggers: Optional[List[Any]] = None
    skipCQR: bool = False

    # diphoton preselection cuts
    min_pt_photon = 25.0
    min_pt_lead_photon = 35.0
    min_mvaid = -0.9
    max_sc_eta = 2.5
    gap_barrel_eta = 1.4442
    gap_endcap_eta = 1.566
    max_hovere = 0.08
    min_full5x5_r9 = 0.8
    max_chad_iso = 20.0
    max_chad_rel_iso = 0.3

    min_full5x5_r9_EB_high_r9 = 0.85
    min_full5x5_r9_EE_high_r9 = 0.9
    min_full5x5_r9_EB_low_r9 = 0.5
    min_full5x5_r9_EE_low_r9 = 0.8
    max_trkSumPtHollowConeDR03_EB_low_r9 = 6.0
    max_trkSumPtHollowConeDR03_EE_low_r9 = 6.0
    max_sieie_EB_low_r9 = 0.015
    max_sieie_EE_low_r9 = 0.035
    max_pho_iso_EB_low_r9 = 4.0
    max_pho_iso_EE_low_r9 = 4.0

    eta_rho_corr = 1.5
    low_eta_rho_corr = 0.16544
    high_eta_rho_corr = 0.13212

    def __post_init__(self) -> None:
        logger.debug(f"Setting up processor with metaconditions: {self.meta}")

        if self.taggers is None:
            self.taggers = []
        self.taggers.sort(key=lambda x: x.priority)  # type: ignore

        self.prefixes = {"pho_lead": "lead", "pho_sublead": "sublead"}

        # build the chained quantile regressions
        self.chained_quantile = None
        if not self.skipCQR:
            try:
                self.chained_quantile: Optional[
                    ChainedQuantileRegression
                ] = ChainedQuantileRegression(**self.meta["PhoIdInputCorrections"])
            except Exception as e:
                warnings.warn(f"Could not instantiate ChainedQuantileRegression: {e}")
        else:
            logger.info("Skipping CQR as required")

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
        self.diphoton_mva = load_bdt(
            os.path.join(
                diphoton_weights_dir, self.meta["flashggDiPhotonMVA"]["weightFile"]
            )
        )

    def photon_preselection(
        self, photons: awkward.Array, events: awkward.Array
    ) -> awkward.Array:
        # hlt-mimicking cuts
        rho = events.fixedGridRhoAll * awkward.ones_like(photons.pt)
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
                photons.trkSumPtHollowConeDR03
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
                photons.trkSumPtHollowConeDR03
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
                | (photons.pfRelIso03_chg < self.max_chad_iso)
                | (photons.pfRelIso03_chg / photons.pt < self.max_chad_rel_iso)
            )
            & (isEB_high_r9 | isEB_low_r9 | isEE_high_r9 | isEE_low_r9)
        ]

    def diphoton_list_to_pandas(self, diphotons: awkward.Array) -> pandas.DataFrame:
        output = pandas.DataFrame()
        for field in awkward.fields(diphotons):
            prefix = self.prefixes.get(field, "")
            if len(prefix) > 0:
                for subfield in awkward.fields(diphotons[field]):
                    logger.debug(f"Adding {prefix}_{subfield}")
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
        subdirs = subdirs or []
        xrd_prefix = "root://"
        pfx_len = len(xrd_prefix)
        xrootd = False
        if xrd_prefix in location:
            try:
                import XRootD
                import XRootD.client

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
            triggers = self.meta["TriggerPaths"][self.trigger_group][self.analysis]
            triggered = functools.reduce(
                operator.or_, (events.HLT[trigger[4:-1]] for trigger in triggers)
            )

        return events[filtered & triggered]

    def process(self, events: awkward.Array) -> Dict[Any, Any]:
        # data or monte carlo?
        self.data_kind = "mc" if "GenPart" in awkward.fields(events) else "data"

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
        # sort photons in each event descending in pt
        # make descending-pt combinations of photons
        photons = photons[awkward.argsort(photons.pt, ascending=False)]
        diphotons = awkward.combinations(photons, 2, fields=["pho_lead", "pho_sublead"])
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

        # set diphotons as part of the event record
        events["diphotons"] = diphotons

        # here we start recording possible coffea accumulators
        # most likely histograms, could be counters, arrays, ...
        histos_etc = {}

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
                events
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
            tags = awkward.from_regular(awkward.unflatten(flat_tags, counts), axis=2)
            winner = awkward.min(tags[tags != 0], axis=2)
            diphotons["best_tag"] = winner

            # lowest priority is most important (ascending sort)
            # leave in order of diphoton pT in case of ties (stable sort)
            sorted = awkward.argsort(diphotons.best_tag, stable=True)
            diphotons = diphotons[sorted]

        diphotons = awkward.firsts(diphotons)

        # annotate diphotons with event information
        diphotons["event"] = events.event
        diphotons["lumi"] = events.luminosityBlock
        diphotons["run"] = events.run

        # drop events without a preselected diphoton candidate
        # drop events without a tag, if there are tags
        if len(self.taggers):
            diphotons = diphotons[
                ~(awkward.is_none(diphotons) | awkward.is_none(diphotons.best_tag))
            ]
        else:
            diphotons = diphotons[~awkward.is_none(diphotons)]

        if self.output_location is not None:
            df = self.diphoton_list_to_pandas(diphotons)
            fname = (
                events.behavior["__events_factory__"]._partition_key.replace("/", "_")
                + ".parquet"
            )
            subdirs = []
            if "dataset" in events.metadata:
                subdirs.append(events.metadata["dataset"])
            self.dump_pandas(df, fname, self.output_location, subdirs)

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
        photons["fixedGridRhoAll"] = events.fixedGridRhoAll * awkward.ones_like(
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
