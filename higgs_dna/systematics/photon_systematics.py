import numpy as np


def photon_pt_scale(pt):
    return (1.0 + np.array([0.05, -0.05], dtype=np.float32)) * pt[:, None]
