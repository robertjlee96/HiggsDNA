from .photon_systematics import photon_pt_scale


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
