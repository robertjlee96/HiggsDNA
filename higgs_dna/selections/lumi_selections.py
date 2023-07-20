from coffea.lumi_tools import LumiMask
import awkward
import os
import logging


def select_lumis(
    self,
    events: awkward.highlevel.Array,
    logger: logging.Logger,
) -> awkward.highlevel.Array:
    goldenJson_dict = {
        "2022": os.path.join(
            os.path.dirname(__file__),
            "../metaconditions/CAF/certification/Collisions22/Cert_Collisions2022_355100_362760_Golden.json",
        ),
        "2023": os.path.join(
            os.path.dirname(__file__),
            "../metaconditions/CAF/certification/Collisions23/Cert_Collisions2023_366442_370092_Golden.json",
        ),
    }
    # Reference
    # https://github.com/CoffeaTeam/coffea/blob/f8a4eb97137e84dd52474d26b8100174da196b57/tests/test_lumi_tools.py
    # https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun3Analysis#Data
    lumimask = LumiMask(goldenJson_dict[self.year])
    logger.debug(
        "Year: {} GoldenJson: {}".format(self.year, goldenJson_dict[self.year])
    )
    return lumimask(events.run, events.luminosityBlock)
