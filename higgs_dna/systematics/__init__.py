from .photon_systematics import photon_pt_scale
from .event_weight_systematics import xSec_DYJets


systematics = {
    "PhotonPtScale": {
        "object": "Photon",
        "args": {
            "kind": "UpDownSystematic",
            "what": "pt",
            "varying_function": photon_pt_scale,
        },
    }
}

# "name": varying_function
weight_systematics = {
    "xSecDYJets": xSec_DYJets
}
