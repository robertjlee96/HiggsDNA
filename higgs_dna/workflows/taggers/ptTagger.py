class ptTagger:

    def __init__(self) -> None:
        pass

    @property
    def name(self) -> str:
        return "ptTagger"

    @property
    def priority(self) -> int:
        return 30

    def GetCategory(self, ievt: int) -> int:

        evt_pt = self.pt[ievt][0]

        if evt_pt < 5:
            cat = 0
        elif evt_pt >= 5 and evt_pt < 10:
            cat = 1
        elif evt_pt >= 10 and evt_pt < 15:
            cat = 2
        elif evt_pt >= 15 and evt_pt < 20:
            cat = 3
        elif evt_pt >= 20 and evt_pt < 25:
            cat = 4
        elif evt_pt >= 25 and evt_pt < 30:
            cat = 5
        elif evt_pt >= 30 and evt_pt < 35:
            cat = 6
        elif evt_pt >= 35 and evt_pt < 45:
            cat = 7
        elif evt_pt >= 45 and evt_pt < 60:
            cat = 8
        elif evt_pt >= 60 and evt_pt < 80:
            cat = 9
        elif evt_pt >= 80 and evt_pt < 100:
            cat = 10
        elif evt_pt >= 100 and evt_pt < 120:
            cat = 11
        elif evt_pt >= 120 and evt_pt < 140:
            cat = 12
        elif evt_pt >= 140 and evt_pt < 170:
            cat = 13
        elif evt_pt >= 170 and evt_pt < 200:
            cat = 14
        elif evt_pt >= 200 and evt_pt < 250:
            cat = 15
        elif evt_pt >= 250 and evt_pt < 350:
            cat = 16
        elif evt_pt >= 350 and evt_pt < 450:
            cat = 17
        else:
            cat = 18

        return cat

    def __call__(self, events: ak.Array, diphotons: ak.Array) -> ak.Array:

        self.pt = events.diphotons.pt

        nDiphotons = ak.num(
            events.diphotons.pt, axis=1
        )

        ievts_by_dipho = ak.flatten(
            ak.Array([nDipho * [evt_i] for evt_i, nDipho in enumerate(nDiphotons)])
        )

        cat_vals = ak.Array(map(self.GetCategory, ievts_by_dipho))
        cats = ak.unflatten(cat_vals, nDiphotons)  # Back to size of events.
        cats_by_diphoEvt = self.priority + cats

        return (cats_by_diphoEvt, {})
