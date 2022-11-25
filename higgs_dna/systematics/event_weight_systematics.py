import numpy as np


def xSec_DYJets(Weight):
    """This is a dummy."""
    Weight.add(
        name="xSec_DYJets",
        weight=np.ones(len(Weight._weight)),  # 1 for "systematics", != 1 for "corrections". This is multiplied by the nominal event weight (Weight._weight)
        weightUp=1.05 * np.ones(len(Weight._weight)),
        weightDown=0.95 * np.ones(len(Weight._weight))
    )
    return Weight
