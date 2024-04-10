import awkward as ak
import numpy as np
from coffea.analysis_tools import PackedSelection
import logging

logger = logging.getLogger(__name__)


class fiducialClassicalTagger:
    def __init__(self) -> None:
        pass

    @property
    def name(self) -> str:
        return "fiducialClassicalTagger"

    # lower priority is better
    # first decimal point is category within tag (in case for untagged)
    @property
    def priority(self) -> int:
        return 15

    def __call__(self, events: ak.Array, diphotons: ak.Array) -> ak.Array:
        """GenIsolatedPhoton is used in this fiducial region definition

        Args:
            events (ak.Array): _description_

        Returns:
            ak.Array:
                priority+0: out of fiducial
                priority+1: in fiducial
        """
        selection = PackedSelection()

        # This tagger determines if an event passes the fiducial selection at gen level
        # Therefore, we need to protect data events from this tagger as they do not have the relevant field
        if 'GenIsolatedPhoton' in events.fields:  # if we are processing MC
            if 'iso' in events.GenPart.fields:  # Check if the branch GenPart_iso is in the nanoAOD
                sel_pho = (events.GenPart.pdgId == 22) & (events.GenPart.status == 1) & (events.GenPart.iso * events.GenPart.pt < 10)
                photons = events.GenPart[sel_pho]
                selection.add("n_iso_photon", ak.num(photons) > 1)
                photons = photons.mask[selection.all("n_iso_photon")]
                photons = photons[ak.argsort(photons.pt, ascending=False)]
                lead_pho = photons[:, 0]
                sublead_pho = photons[:, 1]
            else:
                selection.add("n_iso_photon", ak.num(events.GenIsolatedPhoton, axis=1) > 1)
                events = events.mask[selection.all("n_iso_photon")]
                lead_pho = events.GenIsolatedPhoton[:, 0]
                sublead_pho = events.GenIsolatedPhoton[:, 1]
            diphoton = (lead_pho) + (sublead_pho)
            selection.add("lead_photon_scaled_pt", lead_pho.pt / diphoton.mass > 1 / 3)
            selection.add("sublead_photon_scaled_pt", sublead_pho.pt / diphoton.mass > 1 / 4)
            selection.add(
                "lead_photon_eta",
                (np.abs(lead_pho.eta) < 1.4442)
                | ((np.abs(lead_pho.eta) < 2.5) & (np.abs(lead_pho.eta) > 1.566)),
            )
            selection.add(
                "sublead_photon_eta",
                (np.abs(sublead_pho.eta) < 1.4442)
                | ((np.abs(sublead_pho.eta) < 2.5) & (np.abs(sublead_pho.eta) > 1.566)),
            )

            sel_fiducial = selection.all(
                "n_iso_photon",
                "lead_photon_scaled_pt",
                "sublead_photon_scaled_pt",
                "lead_photon_eta",
                "sublead_photon_eta",
            )
            logger.debug(
                f"""_summary_
                before: \t {len(events)}
                ---
                n_iso_photon: \t {ak.sum(selection.all("n_iso_photon"))}
                lead_photon_scaled_pt: \t {ak.sum(selection.all("lead_photon_scaled_pt"))}
                sublead_photon_scaled_pt: \t {ak.sum(selection.all("sublead_photon_scaled_pt"))}
                lead_photon_eta: \t {ak.sum(selection.all("lead_photon_eta"))}
                sublead_photon_eta: \t {ak.sum(selection.all("sublead_photon_eta"))}
                ---
                total in fiducial: \t {ak.sum(sel_fiducial)}
                -end-
                """
            )
            # priority+1: in fiducial
            # priority+0: out of fiducial
            return (
                self.priority + sel_fiducial * 1,
                {},
            )
        else:  # we are processing data
            logger.debug("Processing data and you want the fiducialClassicalTagger. Thus, code of -99 is applied for all events in data.")
            return (
                ak.ones_like(events.event) * -99,  # -99 is our code for data here
                {},
            )
