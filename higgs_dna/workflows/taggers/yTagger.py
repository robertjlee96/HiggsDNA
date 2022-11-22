import awkward as ak
import numpy as np


class yTagger:

    def __init__(self) -> None:
        pass

    @property
    def name(self) -> str:
        return "yTagger"

    @property
    def priority(self) -> int:
        return 20

    def GetCategory(self, ievt: int) -> int:

        evt_pz = self.pz[ievt][0]
        evt_energy = self.energy[ievt][0]

        evt_y = np.arctanh(evt_pz / evt_energy)

        if abs(evt_y) < 0.1:
            cat = 0
        elif (abs(evt_y) >= 0.1) and (abs(evt_y) < 0.2):
            cat = 1
        elif (abs(evt_y) >= 0.2) and (abs(evt_y) < 0.3):
            cat = 2
        elif (abs(evt_y) >= 0.3) and (abs(evt_y) < 0.45):
            cat = 3
        elif (abs(evt_y) >= 0.45) and (abs(evt_y) < 0.6):
            cat = 4
        elif (abs(evt_y) >= 0.6) and (abs(evt_y) < 0.75):
            cat = 5
        elif (abs(evt_y) >= 0.75) and (abs(evt_y) < 0.90):
            cat = 6
        elif (abs(evt_y) >= 0.90) and (abs(evt_y) < 2.5):
            cat = 7

        return cat

    def __call__(self, events: ak.Array, diphotons: ak.Array) -> ak.Array:

        self.pz = events.diphotons.pz

        self.energy = events.diphotons.energy

        nDiphotons = ak.num(
            events.diphotons.pz, axis=1
        )

        ievts_by_dipho = ak.flatten(
            ak.Array([nDipho * [evt_i] for evt_i, nDipho in enumerate(nDiphotons)])
        )

        cat_vals = ak.Array(map(self.GetCategory, ievts_by_dipho))
        cats = ak.unflatten(cat_vals, nDiphotons)  # Back to size of events.
        cats_by_diphoEvt = self.priority + cats

        return (cats_by_diphoEvt, {})
