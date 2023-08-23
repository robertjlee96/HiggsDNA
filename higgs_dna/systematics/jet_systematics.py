from coffea.jetmet_tools import JetCorrectionUncertainty

# from coffea.jetmet_tools import FactorizedJetCorrector
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory
from coffea.lookup_tools import extractor
import awkward as ak
import numpy as np
import os

# from coffea.nanoevents import NanoEventsFactory, NanoAODSchema


def JERC_jet(pt, events, year="2022postEE", is_correction=True):
    # n_jets = ak.num(events.Jet.pt)
    ext = extractor()
    ext.add_weight_sets(
        [
            "* * {}".format(
                os.path.join(
                    os.path.dirname(__file__),
                    "data/Winter22Run3_V2_MC_L1FastJet_AK4PFPuppi.txt",
                )
            ),
            "* * {}".format(
                os.path.join(
                    os.path.dirname(__file__),
                    "data/Winter22Run3_V2_MC_L2L3Residual_AK4PFPuppi.txt",
                )
            ),
            "* * {}".format(
                os.path.join(
                    os.path.dirname(__file__),
                    "data/Winter22Run3_V2_MC_L2Relative_AK4PFPuppi.txt",
                )
            ),
            "* * {}".format(
                os.path.join(
                    os.path.dirname(__file__),
                    "data/Winter22Run3_V2_MC_L2Residual_AK4PFPuppi.txt",
                )
            ),
            "* * {}".format(
                os.path.join(
                    os.path.dirname(__file__),
                    "data/Winter22Run3_V2_MC_L3Absolute_AK4PFPuppi.txt",
                )
            ),
            "* * {}".format(
                os.path.join(
                    os.path.dirname(__file__),
                    "data/Winter22Run3_V2_MC_Uncertainty_AK4PFPuppi.junc.txt",
                )
            ),
        ]
    )
    # If there are new JEC files, you can always replace them!
    ext.finalize()

    jec_stack_names = [
        "Winter22Run3_V2_MC_L1FastJet_AK4PFPuppi",
        "Winter22Run3_V2_MC_L2L3Residual_AK4PFPuppi",
        "Winter22Run3_V2_MC_L2Relative_AK4PFPuppi",
        "Winter22Run3_V2_MC_L2Residual_AK4PFPuppi",
        "Winter22Run3_V2_MC_L3Absolute_AK4PFPuppi",
        "Winter22Run3_V2_MC_Uncertainty_AK4PFPuppi",
    ]

    evaluator = ext.make_evaluator()

    jec_inputs = {name: evaluator[name] for name in jec_stack_names}
    jec_stack = JECStack(jec_inputs)

    name_map = jec_stack.blank_name_map
    name_map["JetPt"] = "pt"
    name_map["JetMass"] = "mass"
    name_map["JetEta"] = "eta"
    name_map["JetA"] = "area"

    jets = events.Jet

    jets["pt_raw"] = (1 - jets["rawFactor"]) * jets["pt"]
    jets["mass_raw"] = (1 - jets["rawFactor"]) * jets["mass"]
    jets["pt_gen"] = ak.values_astype(ak.fill_none(jets.matched_gen.pt, 0), np.float32)
    jets["rho"] = ak.broadcast_arrays(events.Rho.fixedGridRhoFastjetAll, jets.pt)[0]
    name_map["ptGenJet"] = "pt_gen"
    name_map["ptRaw"] = "pt_raw"
    name_map["massRaw"] = "mass_raw"
    name_map["Rho"] = "rho"

    events_cache = events.caches[0]
    # corrector = FactorizedJetCorrector(
    #     Winter22Run3_V2_MC_L2Relative_AK4PFPuppi=evaluator['Winter22Run3_V2_MC_L2Relative_AK4PFPuppi'],
    # )
    uncertainties = JetCorrectionUncertainty(
        Winter22Run3_V2_MC_Uncertainty_AK4PFPuppi=evaluator[
            "Winter22Run3_V2_MC_Uncertainty_AK4PFPuppi"
        ]
    )

    jet_factory = CorrectedJetsFactory(name_map, jec_stack)

    if is_correction:
        corrected_jets = jet_factory.build(jets, lazy_cache=events_cache)
        events.Jet = corrected_jets

        return events

    else:
        sys = list(uncertainties.getUncertainty(JetEta=jets.eta, JetPt=jets.pt))
        level, corrs = sys[0]
        uncertainty = ak.flatten(corrs)

        return uncertainty * (ak.flatten(jets.pt)[:, None])
