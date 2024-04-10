import awkward as ak
import numpy as np
# from higgs_dna.selections.object_selections import delta_r_mask
import logging

logger = logging.getLogger(__name__)


def get_fiducial_flag(events: ak.Array, flavour: str = "Geometric") -> ak.Array:
    """
    Calculate the fiducial flag for events based on photon kinematics and geometric criteria at particle level.

    The function processes the events, identifying those that meet
    specific criteria based on the properties of the leading and subleading photons.
    The fiducial flag is determined based on transverse momentum (pt), mass, and pseudorapidity (eta)
    of the photon pairs, applying either 'Geometric' or 'Classical' selection criteria.

    Parameters:
    - events (ak.Array): An Awkward Array containing event data with GenIsolatedPhoton fields.
    - flavour (str, optional): The selection criterion to apply. Defaults to "Geometric".
      Can be "Geometric" for geometric mean based selection (https://arxiv.org/abs/2106.08329) or "Classical" for classical CMS pt/mass scaled cuts.

    Returns:
    - ak.Array: An Awkward Array of boolean flags, where True indicates an event meets the fiducial
      selection criteria.

    Note:
    - The function pads GenIsolatedPhoton fields to ensure at least two photons are present per event,
      filling missing values with None.
    - If the GenPart_iso branch is not included in the NanoAOD, the GenIsolatedPhotons collection is used
    """
    if 'iso' in events.GenPart.fields:
        sel_pho = (events.GenPart.pdgId == 22) & (events.GenPart.status == 1) & (events.GenPart.iso * events.GenPart.pt < 10)
        photons = events.GenPart[sel_pho]
        photons = photons[ak.argsort(photons.pt, ascending=False)]
        GenIsolatedPhotons = ak.pad_none(photons, 2)
    else:
        # Extract and pad the gen isolated photons
        GenIsolatedPhotons = events.GenIsolatedPhoton
        GenIsolatedPhotons = ak.pad_none(GenIsolatedPhotons, 2)

    # Separate leading and subleading photons
    lead_pho = GenIsolatedPhotons[:, 0]
    sublead_pho = GenIsolatedPhotons[:, 1]

    # Calculate diphoton system four vector
    diphoton = lead_pho + sublead_pho

    # Apply selection criteria based on the specified flavour
    if flavour == 'Geometric':
        # Geometric mean of pt criterion
        lead_mask = np.sqrt(lead_pho.pt * sublead_pho.pt) / diphoton.mass > 1 / 3
    elif flavour == 'Classical':
        # Classical pt/mass ratio criterion
        lead_mask = lead_pho.pt / diphoton.mass > 1 / 3

    # Subleading photon criterion always the same
    sublead_mask = sublead_pho.pt / diphoton.mass > 1 / 4

    # Pseudorapidity criteria for leading and subleading photons
    # Within tracker acceptance and remove the gap region
    # Note: Based on classical eta, not SC eta here
    lead_eta_mask = (np.abs(lead_pho.eta) < 1.4442) | ((np.abs(lead_pho.eta) < 2.5) & (np.abs(lead_pho.eta) > 1.566))
    sublead_eta_mask = (np.abs(sublead_pho.eta) < 1.4442) | ((np.abs(sublead_pho.eta) < 2.5) & (np.abs(sublead_pho.eta) > 1.566))

    # Combine all selection masks to form the fiducial flag
    fiducial_flag = lead_mask & sublead_mask & lead_eta_mask & sublead_eta_mask
    # Fill None values with False
    # Note: These values result from the padding
    # Only occurs for events that did not have two GenIsolatedPhoton origin
    fiducial_flag = ak.fill_none(fiducial_flag, False)

    return fiducial_flag


def get_NGenJets(events: ak.Array, pt_cut, eta_cut) -> ak.Array:
    GenJets = events.GenJet
    # GenParts = events.GenPart
    # GenPart_is_from_Higgs = events.GenPart[GenParts.genPartIdxMother].pdgId == 25
    # events.GenPart[events.GenPart.pdgId==22 & ]

    # dr_dipho_cut = delta_r_mask(jets, diphotons, self.jet_dipho_min_dr)
    GenJets = GenJets[GenJets.pt > pt_cut]
    GenJets = GenJets[np.abs(GenJets.eta) < eta_cut]
    NGenJets = ak.num(GenJets)
    return NGenJets
