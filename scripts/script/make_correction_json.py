#!/usr/bin/env python
import correctionlib
import correctionlib.schemav2 as cs
import rich
from optparse import OptionParser

# This script is meant to create json files containing the correction and systematic variations as were used in flashgg
# the files are created using correctionlib, copying the binning, corrections and errors from flashgg

usage = "Usage: python %prog [options]"
parser = OptionParser(usage=usage)
parser.add_option(
    "--read",
    dest="read",
    default="",
    help="file to read from, if different from empty it will skip the creation step",
)
parser.add_option(
    "--do_fnuf",
    dest="do_fnuf",
    action="store_true",
    default=False,
    help="create FNUF correction json",
)
parser.add_option(
    "--do_showsh",
    dest="do_showsh",
    action="store_true",
    default=False,
    help="create shower shape correction json",
)
parser.add_option(
    "--do_idmva",
    dest="do_IDMVA",
    action="store_true",
    default=False,
    help="create PhotonID MVA correction json",
)
(opt, args) = parser.parse_args()


def multibinning(inputs_: list, edges_: list, content_, flow_: str):
    return cs.MultiBinning(
        nodetype="multibinning",
        inputs=inputs_,
        edges=edges_,
        content=content_,
        flow=flow_,
    )


# FNUF variations fromflashgg: https://github.com/cms-analysis/flashgg/blob/58c243d9d1f794d7dca8a94fdcd390aed91cb49c/Systematics/python/flashggDiPhotonSystematics2017_Legacy_cfi.py#L73-L82
# 4 bins (2 in eta 2 in r9)
# the code creates the CorrectionSet using correctionlib, write it to a .json and reads it back
# in principle one can dump more than one Correction in the same .json adding it to the CorrectionSet
if opt.do_fnuf:
    inputs_ = ["eta", "r9"]
    edges_ = [[0.0, 1.5, 6.0], [0.0, 0.94, 999.0]]
    content_ = {
        "nominal": [1.0, 1.0, 1.0, 1.0],
        "up": [1.0007, 1.0022, 1.00005, 1.00251],
        "down": [0.9993, 0.9978, 0.99995, 0.99749],
    }
    flow_ = "clamp"

    FNUF = cs.Correction(
        name="FNUF",
        version=1,
        inputs=[
            cs.Variable(
                name="systematic", type="string", description="Systematic variation"
            ),
            cs.Variable(name="eta", type="real", description="Photon eta"),
            cs.Variable(
                name="r9",
                type="real",
                description="Photon full 5x5 R9, ratio E3x3/ERAW, where E3x3 is the energy sum of the 3 by 3 crystals surrounding the supercluster seed crystal and ERAW is the raw energy sum of the supercluster",
            ),
        ],
        output=cs.Variable(
            name="Ecorr",
            type="real",
            description="Multiplicative correction to photon energy",
        ),
        data=cs.Category(
            nodetype="category",
            input="systematic",
            content=[
                {
                    "key": "nominal",
                    "value": multibinning(inputs_, edges_, content_["nominal"], flow_),
                },
                {
                    "key": "up",
                    "value": multibinning(inputs_, edges_, content_["up"], flow_),
                },
                {
                    "key": "down",
                    "value": multibinning(inputs_, edges_, content_["down"], flow_),
                },
            ],
        ),
    )

    rich.print(FNUF)

    # test that everything is fine with some mockup values
    etas = [1.0, 2.0, 1.0, 2.0]
    rs = [0.9, 0.99, 0.99, 0.9]

    print("-" * 120)
    print("nominal --->","eta:",etas,"r9:",rs,"result:",FNUF.to_evaluator().evaluate("nominal", etas, rs))
    print("-" * 120)
    print("up      --->","eta:",etas,"r9:",rs,"result:",FNUF.to_evaluator().evaluate("up", etas, rs))
    print("-" * 120)
    print("down    --->","eta:",etas,"r9:",rs,"result:",FNUF.to_evaluator().evaluate("down", etas, rs))
    print("-" * 120)
    print()

    cset = cs.CorrectionSet(
        schema_version=2,
        description="FNUF",
        corrections=[
            FNUF,
        ],
    )

    # if we're not just checking an existing json we create a new one
    if opt.read == "":
        with open("FNUF.json", "w") as fout:
            fout.write(cset.json(exclude_unset=True))

        import gzip

        with gzip.open("FNUF.json.gz", "wt") as fout:
            fout.write(cset.json(exclude_unset=True))

    print("Reading back...")

    file = opt.read
    ceval = correctionlib.CorrectionSet.from_file(file)
    for corr in ceval.values():
        print(f"Correction {corr.name} has {len(corr.inputs)} inputs")
        for ix in corr.inputs:
            print(f"   Input {ix.name} ({ix.type}): {ix.description}")

    print("-" * 120)
    print("nominal --->","eta:",etas,"r9:",rs,"result:",ceval["FNUF"].evaluate("nominal", etas, rs))
    print("-" * 120)
    print("up      --->","eta:",etas,"r9:",rs,"result:",ceval["FNUF"].evaluate("up", etas, rs))
    print("-" * 120)
    print("down    --->","eta:",etas,"r9:",rs,"result:",ceval["FNUF"].evaluate("down", etas, rs))
    print("-" * 120)

# Shower Shape variations fromflashgg: https://github.com/cms-analysis/flashgg/blob/58c243d9d1f794d7dca8a94fdcd390aed91cb49c/Systematics/python/flashggDiPhotonSystematics2017_Legacy_cfi.py#L73-L82
# 8 bins (4 in eta 2 in r9)
# Flashgg correction is from Martina: https://indico.cern.ch/event/628676/contributions/2546615/attachments/1440085/2216643/20170405_martina_regrEchecksUpdate.pdf

if opt.do_showsh:
    inputs_ = ["eta", "r9"]
    edges_ = [[0.0, 1.0, 1.5, 2.0, 6.0], [0.0, 0.94, 999.0]]
    content_ = {
        "nominal": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        "up": [0.9999, 0.9994, 1.0002, 0.9989, 1.0015, 1.0002, 1.0004, 1.0003],
        "down": [1.0001, 1.0006, 0.9998, 1.0011, 0.9985, 0.9998, 0.9996, 0.9997],
    }
    flow_ = "clamp"

    SS = cs.Correction(
        name="ShowerShape",
        version=1,
        inputs=[
            cs.Variable(
                name="systematic", type="string", description="Systematic variation"
            ),
            cs.Variable(name="eta", type="real", description="Photon eta"),
            cs.Variable(
                name="r9",
                type="real",
                description="Photon full 5x5 R9, ratio E3x3/ERAW, where E3x3 is the energy sum of the 3 by 3 crystals surrounding the supercluster seed crystal and ERAW is the raw energy sum of the supercluster",
            ),
        ],
        output=cs.Variable(
            name="Ecorr",
            type="real",
            description="Multiplicative correction to photon energy and pt",
        ),
        data=cs.Category(
            nodetype="category",
            input="systematic",
            content=[
                {
                    "key": "nominal",
                    "value": multibinning(inputs_, edges_, content_["nominal"], flow_),
                },
                {
                    "key": "up",
                    "value": multibinning(inputs_, edges_, content_["up"], flow_),
                },
                {
                    "key": "down",
                    "value": multibinning(inputs_, edges_, content_["down"], flow_),
                },
            ],
        ),
    )

    rich.print(SS)

    etas = [0.8,1.2,1.7,2.5,0.8,1.2,1.7,2.0]
    rs = [0.9, 0.9, 0.9, 0.9, 0.99, 0.99, 0.99, 0.99]

    print("-" * 120)
    print("nominal --->","eta:",etas,"r9:",rs,"result:",SS.to_evaluator().evaluate("nominal", etas, rs))
    print("-" * 120)
    print("up      --->","up      --->","eta:",etas,"r9:",rs,"result:",SS.to_evaluator().evaluate("up", etas, rs))
    print("-" * 120)
    print("down    --->","eta:",etas,"r9:",rs,"result:",SS.to_evaluator().evaluate("down", etas, rs))
    print("-" * 120)
    print()

    cset = cs.CorrectionSet(
        schema_version=2,
        description="ShowerShape",
        corrections=[
            SS,
        ],
    )

    if opt.read == "":
        with open("ShowerShape.json", "w") as fout:
            fout.write(cset.json(exclude_unset=True))

        import gzip

        with gzip.open("ShowerShape.json.gz", "wt") as fout:
            fout.write(cset.json(exclude_unset=True))

    if opt.read != "":
        print("Reading back...")
        file = opt.read
        ceval = correctionlib.CorrectionSet.from_file(file)
        for corr in ceval.values():
            print(f"Correction {corr.name} has {len(corr.inputs)} inputs")
            for ix in corr.inputs:
                print(f"   Input {ix.name} ({ix.type}): {ix.description}")

        print("-" * 120)
        print("nominal --->","eta:",etas,"r9:",rs,"result:",ceval["ShowerShape"].evaluate("nominal", etas, rs))
        print("-" * 120)
        print("up      --->","eta:",etas,"r9:",rs,"result:",ceval["ShowerShape"].evaluate("up", etas, rs))
        print("-" * 120)
        print("down    --->","eta:",etas,"r9:",rs,"result:",ceval["ShowerShape"].evaluate("down", etas, rs))
        print("-" * 120)

# PhotonID MVA variations fromflashgg: https://github.com/cms-analysis/flashgg/blob/ae6563050722bd168545eac2b860ef56cdda7be4/Systematics/python/flashggDiPhotonSystematics2017_Legacy_cfi.py#L280-286
# 8 bins (2 in eta 3 in r9)
# LooseIDMVA [-0.9] SF and uncertainty for UL2017. Dt: 17/11/2020
# link to the presentation: https://indico.cern.ch/event/963617/contributions/4103623/attachments/2141570/3608645/Zee_Validation_UL2017_Update_09112020_Prasant.pdf

if opt.do_IDMVA:
    inputs_ = ["eta", "r9"]
    edges_ = [[0.0, 1.5, 6.0], [0.0, 0.85, 0.9, 999.0]]
    content_ = {
        "nominal": [1.0021, 1.0001, 1.0001, 1.0061, 1.0061, 1.0016],
        "up": [1.0035, 1.0039, 1.0039, 1.0075, 1.0075, 1.0039],
        "down": [1.0007, 0.9963, 0.9963, 1.0047, 1.0047, 0.9993],
    }
    flow_ = "clamp"

    MVA = cs.Correction(
        name="LooseMvaSF",
        version=1,
        inputs=[
            cs.Variable(
                name="systematic", type="string", description="Systematic variation"
            ),
            cs.Variable(name="eta", type="real", description="Photon eta"),
            cs.Variable(
                name="r9",
                type="real",
                description="Photon full 5x5 R9, ratio E3x3/ERAW, where E3x3 is the energy sum of the 3 by 3 crystals surrounding the supercluster seed crystal and ERAW is the raw energy sum of the supercluster",
            ),
        ],
        output=cs.Variable(
            name="Ecorr",
            type="real",
            description="Multiplicative correction to photon energy and pt",
        ),
        data=cs.Category(
            nodetype="category",
            input="systematic",
            content=[
                {
                    "key": "nominal",
                    "value": multibinning(inputs_, edges_, content_["nominal"], flow_),
                },
                {
                    "key": "up",
                    "value": multibinning(inputs_, edges_, content_["up"], flow_),
                },
                {
                    "key": "down",
                    "value": multibinning(inputs_, edges_, content_["down"], flow_),
                },
            ],
        ),
    )

    rich.print(MVA)

    etas = [0.8, 0.8, 0.8, 2.5, 2.5, 2.5]
    rs = [0.6, 0.87, 0.95, 0.6, 0.87, 0.99]

    print("-" * 120)
    print("nominal --->","eta:",etas,"r9:",rs,"result:",MVA.to_evaluator().evaluate("nominal", etas, rs))
    print("-" * 120)
    print("up      --->","eta:",etas,"r9:",rs,"result:",MVA.to_evaluator().evaluate("up", etas, rs))
    print("-" * 120)
    print("down    --->","eta:",etas,"r9:",rs,"result:",MVA.to_evaluator().evaluate("down", etas, rs))
    print("-" * 120)
    print()

    cset = cs.CorrectionSet(
        schema_version=2,
        description="LooseMvaSF",
        corrections=[
            MVA,
        ],
    )

    if opt.read == "":
        with open("LooseMvaSF.json", "w") as fout:
            fout.write(cset.json(exclude_unset=True))

        import gzip

        with gzip.open("LooseMvaSF.json.gz", "wt") as fout:
            fout.write(cset.json(exclude_unset=True))

    if opt.read != "":
        print("Reading back...")
        file = opt.read
        ceval = correctionlib.CorrectionSet.from_file(file)
        for corr in ceval.values():
            print(f"Correction {corr.name} has {len(corr.inputs)} inputs")
            for ix in corr.inputs:
                print(f"   Input {ix.name} ({ix.type}): {ix.description}")

        print("-" * 120)
        print("nominal --->","eta:",etas,"r9:",rs,"result:",ceval["LooseMvaSF"].evaluate("nominal", etas, rs))
        print("-" * 120)
        print("up      --->","eta:",etas,"r9:",rs,"result:",ceval["LooseMvaSF"].evaluate("up", etas, rs))
        print("-" * 120)
        print("down    --->","eta:",etas,"r9:",rs,"result:",ceval["LooseMvaSF"].evaluate("down", etas, rs))
        print("-" * 120)
