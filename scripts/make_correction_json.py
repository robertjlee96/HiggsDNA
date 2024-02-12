#!/usr/bin/env python
import correctionlib
import correctionlib.schemav2 as cs
import rich
from optparse import OptionParser
import pydantic

print(pydantic.__version__)
print(correctionlib.__version__)

# This script is meant to create json files containing the correction and systematic variations as were used in flashgg
# the files are created using correctionlib, copying the binning, corrections and errors from flashgg

# The script works with the most recent version of correctionlib, to be compatible with older releases one should 
# change the istances of "cset.model_dump_json" to "cset.json" to not stumble into errors

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
parser.add_option(
    "--do_Mat",
    dest="do_Material",
    action="store_true",
    default=False,
    help="create Material correction json",
)
parser.add_option(
    "--do_eVeto",
    dest="do_eVeto",
    action="store_true",
    default=False,
    help="create Electron Veto correction json",
)
parser.add_option(
    "--do_presel",
    dest="do_presel",
    action="store_true",
    default=False,
    help="create Preselection Scale Factor json",
)
parser.add_option(
    "--do_trigger",
    dest="do_trigger",
    action="store_true",
    default=False,
    help="create Trigger Scale Factor jsons",
)
parser.add_option(
    "--do_trigger_lead",
    dest="do_trigger_lead",
    action="store_true",
    default=False,
    help="create lead photon Trigger Scale Factor json",
)
parser.add_option(
    "--do_trigger_sublead",
    dest="do_trigger_sublead",
    action="store_true",
    default=False,
    help="create sublead photon Trigger Scale Factor json",
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

# Material  variations fromflashgg: https://github.com/cms-analysis/flashgg/blob/ae6563050722bd168545eac2b860ef56cdda7be4/Systematics/python/flashggDiPhotonSystematics2017_Legacy_cfi.py#L280-286
# 6 bins (3 in SCeta 2 in r9)

if opt.do_Material:
    inputs_ = ["SCEta", "r9"]
    edges_ = [[0.0, 1., 1.5, 999.0], [0.0, 0.94, 999.0]]
    content_ = {
        "nominal": [1., 1., 1., 1., 1., 1.],
        "up": [1.000455, 1.000233, 1.002089, 1.002089, 1.001090, 1.002377],
        "down": [0.999545, 0.999767, 0.9979911, 0.9979911, 0.99891, 0.997623],
    }
    flow_ = "clamp"

    MAT = cs.Correction(
        name="Material",
        description="Material correction",
        generic_formulas=None,
        version=1,
        inputs=[
            cs.Variable(
                name="systematic", type="string", description="Systematic variation"
            ),
            cs.Variable(name="SCEta", type="real", description="Photon Super Cluster eta"),
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
                cs.CategoryItem(
                    key="up",
                    value=multibinning(inputs_, edges_, content_["up"], flow_)
                ),
                cs.CategoryItem(
                    key="down",
                    value=multibinning(inputs_, edges_, content_["down"], flow_)
                ),
            ],
            default=multibinning(inputs_, edges_, content_["nominal"], flow_),
        ),
    )

    rich.print(MAT)

    etas = [0.8, 0.8, 1.2, 1.2, 2.5, 2.5]
    rs = [0.6, 0.95, 0.6, 0.95, 0.6, 0.99]

    print(MAT.to_evaluator())
    print("-" * 120)
    print("nominal --->","eta:",etas,"r9:",rs,"result:",MAT.to_evaluator().evaluate("nominal", etas, rs))
    print("-" * 120)
    print("up      --->","eta:",etas,"r9:",rs,"result:",MAT.to_evaluator().evaluate("up", etas, rs))
    print("-" * 120)
    print("down    --->","eta:",etas,"r9:",rs,"result:",MAT.to_evaluator().evaluate("down", etas, rs))
    print("-" * 120)
    print()

    cset = cs.CorrectionSet(
        schema_version=2,
        description="Material",
        compound_corrections=None,
        corrections=[
            MAT,
        ],
    )

    if opt.read == "":
        with open("Material.json", "w") as fout:
            fout.write(cset.model_dump_json(exclude_unset=True))

        import gzip

        with gzip.open("Material.json.gz", "wt") as fout:
            fout.write(cset.model_dump_json(exclude_unset=True))

    if opt.read != "":
        print("Reading back...")
        file = opt.read
        ceval = correctionlib.CorrectionSet.from_file(file)
        for corr in ceval.values():
            print(f"Correction {corr.name} has {len(corr.inputs)} inputs")
            for ix in corr.inputs:
                print(f"   Input {ix.name} ({ix.type}): {ix.description}")

        print("-" * 120)
        print("nominal --->","eta:",etas,"r9:",rs,"result:",ceval["Material"].evaluate("nominal", etas, rs))
        print("-" * 120)
        print("up      --->","eta:",etas,"r9:",rs,"result:",ceval["Material"].evaluate("up", etas, rs))
        print("-" * 120)
        print("down    --->","eta:",etas,"r9:",rs,"result:",ceval["Material"].evaluate("down", etas, rs))
        print("-" * 120)


# elecronVetoSF  variations fromflashgg: https://github.com/cms-analysis/flashgg/blob/4edea8897e2a4b0518dca76ba6c9909c20c40ae7/Systematics/python/flashggDiPhotonSystematics2017_Legacy_cfi.py#L27
# 6 bins (2 in SCeta 3 in r9)
# link to the presentation: JTao: slide 13 of the updated UL2017 results, https://indico.cern.ch/event/961164/contributions/4089584/attachments/2135019/3596299/Zmmg_UL2017%20With%20CorrMC_Hgg%20%2802.11.2020%29.pdf, presented by Aamir on 2nd Nov. 2020

if opt.do_eVeto:
    inputs_ = ["SCEta", "r9"]
    edges_ = [[0.0, 1.5, 6.0], [0.0, 0.85, 0.9, 999.0]]
    content_ = {
        "nominal": [0.9838, 0.9913, 0.9913, 0.9777, 0.9777, 0.9784],
        "up": [0.9862, 0.9922, 0.9922, 0.9957, 0.9957, 0.981],
        "down": [0.9814, 0.9904, 0.9904, 0.9597, 0.9597, 0.9758],
    }
    flow_ = "clamp"

    EVETO = cs.Correction(
        name="ElectronVetoSF",
        description="Electron Veto Scale Factor",
        generic_formulas=None,
        version=1,
        inputs=[
            cs.Variable(
                name="systematic", type="string", description="Systematic variation"
            ),
            cs.Variable(name="SCEta", type="real", description="Photon Super Cluster eta absolute value"),
            cs.Variable(
                name="r9",
                type="real",
                description="Photon full 5x5 R9, ratio E3x3/ERAW, where E3x3 is the energy sum of the 3 by 3 crystals surrounding the supercluster seed crystal and ERAW is the raw energy sum of the supercluster",
            ),
        ],
        output=cs.Variable(
            name="Ecorr",
            type="real",
            description="Multiplicative correction to event weight (per-photon)",
        ),
        data=cs.Category(
            nodetype="category",
            input="systematic",
            content=[
                cs.CategoryItem(
                    key="up",
                    value=multibinning(inputs_, edges_, content_["up"], flow_)
                ),
                cs.CategoryItem(
                    key="down",
                    value=multibinning(inputs_, edges_, content_["down"], flow_)
                ),
            ],
            default=multibinning(inputs_, edges_, content_["nominal"], flow_),
        ),
    )

    rich.print(EVETO)

    etas = [0.8, 0.8, 0.8, 1.7, 1.7, 1.7]
    rs = [0.6, 0.85, 0.95, 0.6, 0.85, 0.99]

    print(EVETO.to_evaluator())
    print("-" * 120)
    print("nominal --->","eta:",etas,"r9:",rs,"result:",EVETO.to_evaluator().evaluate("nominal", etas, rs))
    print("-" * 120)
    print("up      --->","eta:",etas,"r9:",rs,"result:",EVETO.to_evaluator().evaluate("up", etas, rs))
    print("-" * 120)
    print("down    --->","eta:",etas,"r9:",rs,"result:",EVETO.to_evaluator().evaluate("down", etas, rs))
    print("-" * 120)
    print()

    cset = cs.CorrectionSet(
        schema_version=2,
        description="Electron Veto SF",
        compound_corrections=None,
        corrections=[
            EVETO,
        ],
    )

    if opt.read == "":
        with open("eVetoSF.json", "w") as fout:
            fout.write(cset.model_dump_json(exclude_unset=True))

        import gzip

        with gzip.open("eVetoSF.json.gz", "wt") as fout:
            fout.write(cset.model_dump_json(exclude_unset=True))

    if opt.read != "":
        print("Reading back...")
        file = opt.read
        ceval = correctionlib.CorrectionSet.from_file(file)
        for corr in ceval.values():
            print(f"Correction {corr.name} has {len(corr.inputs)} inputs")
            for ix in corr.inputs:
                print(f"   Input {ix.name} ({ix.type}): {ix.description}")

        print("-" * 120)
        print("nominal --->","eta:",etas,"r9:",rs,"result:",ceval["ElectronVetoSF"].evaluate("nominal", etas, rs))
        print("-" * 120)
        print("up      --->","eta:",etas,"r9:",rs,"result:",ceval["ElectronVetoSF"].evaluate("up", etas, rs))
        print("-" * 120)
        print("down    --->","eta:",etas,"r9:",rs,"result:",ceval["ElectronVetoSF"].evaluate("down", etas, rs))
        print("-" * 120)


# Presel SF and uncertainty for UL2017. Dt:17/11/2020 fromflashgg: https://github.com/cms-analysis/flashgg/blob/4edea8897e2a4b0518dca76ba6c9909c20c40ae7/Systematics/python/flashggDiPhotonSystematics2017_Legacy_cfi.py#L27
# 6 bins (2 in SCeta 3 in r9)
# link to the presentation: https://indico.cern.ch/event/963617/contributions/4103623/attachments/2141570/3608645/Zee_Validation_UL2017_Update_09112020_Prasant.pdf

if opt.do_presel:
    inputs_ = ["SCEta", "r9"]
    edges_ = [[0.0, 1.5, 6.0], [0.0, 0.85, 0.9, 999.0]]
    content_ = {
        "nominal": [0.9961, 0.9981, 0.9981, 1.0054, 1.0054, 1.0061],
        "up": [1.0268, 1.0038, 1.0038, 1.0183, 1.0183, 1.0079],
        "down": [0.9654, 0.9924, 0.9924, 0.9925, 0.9925, 1.0043],
    }
    flow_ = "clamp"

    PRES = cs.Correction(
        name="PreselSF",
        description="Preselection Scale Factor",
        generic_formulas=None,
        version=1,
        inputs=[
            cs.Variable(
                name="systematic", type="string", description="Systematic variation"
            ),
            cs.Variable(name="SCEta", type="real", description="Photon Super Cluster eta absolute value"),
            cs.Variable(
                name="r9",
                type="real",
                description="Photon full 5x5 R9, ratio E3x3/ERAW, where E3x3 is the energy sum of the 3 by 3 crystals surrounding the supercluster seed crystal and ERAW is the raw energy sum of the supercluster",
            ),
        ],
        output=cs.Variable(
            name="Ecorr",
            type="real",
            description="Multiplicative correction to event weight (per-photon)",
        ),
        data=cs.Category(
            nodetype="category",
            input="systematic",
            content=[
                cs.CategoryItem(
                    key="up",
                    value=multibinning(inputs_, edges_, content_["up"], flow_)
                ),
                cs.CategoryItem(
                    key="down",
                    value=multibinning(inputs_, edges_, content_["down"], flow_)
                ),
            ],
            default=multibinning(inputs_, edges_, content_["nominal"], flow_),
        ),
    )

    rich.print(PRES)

    etas = [0.8, 0.8, 0.8, 1.7, 1.7, 1.7]
    rs = [0.6, 0.85, 0.95, 0.6, 0.85, 0.99]

    print(PRES.to_evaluator())
    print("-" * 120)
    print("nominal --->","eta:",etas,"r9:",rs,"result:",PRES.to_evaluator().evaluate("nominal", etas, rs))
    print("-" * 120)
    print("up      --->","eta:",etas,"r9:",rs,"result:",PRES.to_evaluator().evaluate("up", etas, rs))
    print("-" * 120)
    print("down    --->","eta:",etas,"r9:",rs,"result:",PRES.to_evaluator().evaluate("down", etas, rs))
    print("-" * 120)
    print()

    cset = cs.CorrectionSet(
        schema_version=2,
        description="Preselection SF",
        compound_corrections=None,
        corrections=[
            PRES,
        ],
    )

    if opt.read == "":
        with open("PreselSF.json", "w") as fout:
            fout.write(cset.model_dump_json(exclude_unset=True))

        import gzip

        with gzip.open("PreselSF.json.gz", "wt") as fout:
            fout.write(cset.model_dump_json(exclude_unset=True))

    if opt.read != "":
        print("Reading back...")
        file = opt.read
        ceval = correctionlib.CorrectionSet.from_file(file)
        for corr in ceval.values():
            print(f"Correction {corr.name} has {len(corr.inputs)} inputs")
            for ix in corr.inputs:
                print(f"   Input {ix.name} ({ix.type}): {ix.description}")

        print("-" * 120)
        print("nominal --->","eta:",etas,"r9:",rs,"result:",ceval["PreselSF"].evaluate("nominal", etas, rs))
        print("-" * 120)
        print("up      --->","eta:",etas,"r9:",rs,"result:",ceval["PreselSF"].evaluate("up", etas, rs))
        print("-" * 120)
        print("down    --->","eta:",etas,"r9:",rs,"result:",ceval["PreselSF"].evaluate("down", etas, rs))
        print("-" * 120)


# trigger SF and uncertainty for UL2017 lead photon. Dt:17/11/2020 fromflashgg: https://github.com/cms-analysis/flashgg/blob/4edea8897e2a4b0518dca76ba6c9909c20c40ae7/Systematics/python/flashggDiPhotonSystematics2017_Legacy_cfi.py#L27
# a lot of bins (81, 3 in SCeta, 3 in r9 and 9 in Pt)
# link to the presentation: https://indico.cern.ch/event/963617/contributions/4103623/attachments/2141570/3608645/Zee_Validation_UL2017_Update_09112020_Prasant.pdf

if opt.do_trigger_lead or opt.do_trigger:
    inputs_ = ["SCEta", "r9", "pt"]
    edges_ = [[0.0, 1.5, 3, 999.0], [0.0, 0.56, 0.85, 0.9, 999.0], [0., 35., 37., 40., 45., 50., 60., 70., 90., 999999.0]]
    content_ = {
        "nominal": [ #  35            37            40            45            50            60            70            90             inf    pt
            0.5962935611, 0.7301764945, 0.7481329080, 0.7715837809, 0.7885084282, 0.8106800401, 0.8403040278, 0.8403394687, 0.9294116662, # 0 < eta < 1.5, R9 < 0.56
            0.8185203301, 0.9430487014, 0.9502315420, 0.9569881714, 0.9607761449, 0.9667723989, 0.9735556426, 0.9779985569, 0.9846204583, # 0 < eta < 1.5, 0.56 < R9 < 0.85
            0.8487571377, 0.9553428958, 0.9624070784, 0.9693638237, 0.9755313942, 0.9777963391, 0.9812161610, 0.9845542680, 0.9902588867, # 0 < eta < 1.5, 0.85 < R9 < 0.90
            0.8487571377, 0.9553428958, 0.9624070784, 0.9693638237, 0.9755313942, 0.9777963391, 0.9812161610, 0.9845542680, 0.9902588867, # 0 < eta < 1.5, R9 > 0.90
            0.6297968504, 0.8049372754, 0.8314952358, 0.8544229767, 0.8875746672, 0.9033407955, 0.9207605401, 0.9410420565, 0.9586907211, # 1.5 < eta < 3, R9 < 0.56
            0.6297968504, 0.8049372754, 0.8314952358, 0.8544229767, 0.8875746672, 0.9033407955, 0.9207605401, 0.9410420565, 0.9586907211, # 1.5 < eta < 3, 0.56 < R9 < 0.85
            0.7816089799, 0.9601546944, 0.9728943976, 0.9787293111, 0.9836865868, 0.9845440645, 0.9863780801, 0.9913050524, 0.9969391106, # 1.5 < eta < 3, 0.85 < R9 < 0.90
            0.7758811194, 0.9592830235, 0.9692856214, 0.9763703079, 0.9814613177, 0.9825431442, 0.9857720941, 0.9904181104, 0.9923572396, # 1.5 < eta < 3, R9 > 0.90
            1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, # eta > 3, R9 < 0.56
            1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, # eta > 3, 0.56 < R9 < 0.85
            1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, # eta > 3, 0.85 < R9 < 0.90
            1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000  # eta > 3, R9 > 0.90
            ],
        "up": [ #  35            37            40            45            50            60            70            90             inf    pt
            0.6005737887, 0.7328812488, 0.7503826099, 0.7742611691, 0.7919288403, 0.8258226645, 0.8589755363, 0.8629941224, 0.9561358245, # 0 < eta < 1.5, R9 < 0.56
            0.8203873280, 0.9441263182, 0.9512360353, 0.9580805517, 0.9617771241, 0.9688333341, 0.9754818208, 0.9794625844, 0.9859530068, # 0 < eta < 1.5, 0.56 < R9 < 0.85
            0.8578682281, 0.9577316838, 0.9653803117, 0.9703944882, 0.9765665213, 0.9799697256, 0.9858818495, 0.9855805535, 0.9913562563, # 0 < eta < 1.5, 0.85 < R9 < 0.90
            0.8578682281, 0.9577316838, 0.9653803117, 0.9703944882, 0.9765665213, 0.9799697256, 0.9858818495, 0.9855805535, 0.9913562563, # 0 < eta < 1.5, R9 > 0.90
            0.6324307687, 0.8076131776, 0.8336690619, 0.8559467838, 0.8904000969, 0.9114569193, 0.9278738793, 0.9485921217, 0.9650381006, # 1.5 < eta < 3, R9 < 0.56
            0.6324307687, 0.8076131776, 0.8336690619, 0.8559467838, 0.8904000969, 0.9114569193, 0.9278738793, 0.9485921217, 0.9650381006, # 1.5 < eta < 3, 0.56 < R9 < 0.85
            0.7838966125, 0.9615399325, 0.9739169322, 0.9798026624, 0.9848531139, 0.9867883788, 0.9890037354, 0.9930929945, 0.9990857728, # 1.5 < eta < 3, 0.85 < R9 < 0.90
            0.7785424680, 0.9615285198, 0.9703116690, 0.9774045805, 0.9825372871, 0.9841343722, 0.9872986285, 0.9914973187, 0.9960040747, # 1.5 < eta < 3, R9 > 0.90
            1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, # eta > 3, R9 < 0.56
            1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, # eta > 3, 0.56 < R9 < 0.85
            1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, # eta > 3, 0.85 < R9 < 0.90
            1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, # eta > 3, R9 > 0.90
        ],
        "down": [ #  35            37            40            45            50            60            70            90             inf    pt
            0.5920133335, 0.7274717402, 0.7458832061, 0.7689063927, 0.7850880161, 0.7955374157, 0.8216325193, 0.8176848150, 0.9026875079, # 0 < eta < 1.5, R9 < 0.56
            0.8166533322, 0.9419710846, 0.9492270487, 0.9558957911, 0.9597751657, 0.9647114637, 0.9716294644, 0.9765345294, 0.9832879098, # 0 < eta < 1.5, 0.56 < R9 < 0.85
            0.8396460473, 0.9529541078, 0.9594338451, 0.9683331592, 0.9744962671, 0.9756229526, 0.9765504725, 0.9835279825, 0.9891615171, # 0 < eta < 1.5, 0.85 < R9 < 0.90
            0.8396460473, 0.9529541078, 0.9594338451, 0.9683331592, 0.9744962671, 0.9756229526, 0.9765504725, 0.9835279825, 0.9891615171, # 0 < eta < 1.5, R9 > 0.90
            0.6271629321, 0.8022613732, 0.8293214097, 0.8528991696, 0.8847492375, 0.8952246717, 0.9136472009, 0.9334919913, 0.9523433416, # 1.5 < eta < 3, R9 < 0.56
            0.6271629321, 0.8022613732, 0.8293214097, 0.8528991696, 0.8847492375, 0.8952246717, 0.9136472009, 0.9334919913, 0.9523433416, # 1.5 < eta < 3, 0.56 < R9 < 0.85
            0.7793213473, 0.9587694563, 0.9718718630, 0.9776559598, 0.9825200597, 0.9822997502, 0.9837524248, 0.9895171103, 0.9947924484, # 1.5 < eta < 3, 0.85 < R9 < 0.90
            0.7732197708, 0.9570375272, 0.9682595738, 0.9753360353, 0.9803853483, 0.9809519162, 0.9842455597, 0.9893389021, 0.9887104045, # 1.5 < eta < 3, R9 > 0.90
            1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, # eta > 3, R9 < 0.56
            1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, # eta > 3, 0.56 < R9 < 0.85
            1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, # eta > 3, 0.85 < R9 < 0.90
            1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, # eta > 3, R9 > 0.90
        ]
    }
    flow_ = "clamp"

    TRIG = cs.Correction(
        name="TriggerSF",
        description="Lead Photon Trigger Scale Factor",
        generic_formulas=None,
        version=1,
        inputs=[
            cs.Variable(
                name="systematic", type="string", description="Systematic variation"
            ),
            cs.Variable(name="SCEta", type="real", description="Photon Super Cluster eta absolute value"),
            cs.Variable(
                name="r9",
                type="real",
                description="Photon full 5x5 R9, ratio E3x3/ERAW, where E3x3 is the energy sum of the 3 by 3 crystals surrounding the supercluster seed crystal and ERAW is the raw energy sum of the supercluster",
            ),
            cs.Variable(
                name="pt",
                type="real",
                description="Photon pt",
            ),
        ],
        output=cs.Variable(
            name="Wcorr",
            type="real",
            description="Multiplicative correction to event weight (per-photon)",
        ),
        data=cs.Category(
            nodetype="category",
            input="systematic",
            content=[
                cs.CategoryItem(
                    key="up",
                    value=multibinning(inputs_, edges_, content_["up"], flow_)
                ),
                cs.CategoryItem(
                    key="down",
                    value=multibinning(inputs_, edges_, content_["down"], flow_)
                ),
            ],
            default=multibinning(inputs_, edges_, content_["nominal"], flow_),
        ),
    )

    rich.print(TRIG)

    etas = [
        0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
        0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
        0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
        0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
        1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7,
        1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7,
        1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7,
        1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7,
        10., 10., 10., 10., 10., 10., 10., 10., 10.,
        10., 10., 10., 10., 10., 10., 10., 10., 10.,
        10., 10., 10., 10., 10., 10., 10., 10., 10.,
        10., 10., 10., 10., 10., 10., 10., 10., 10.
        ]
    rs = [
        0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
        0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
        0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88,
        0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99,
        0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
        0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
        0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88,
        0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99,
        0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
        0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
        0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88,
        0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99
        ]
    pts = [
        15., 36., 38., 42., 47., 55., 65., 75., 200.,
        15., 36., 38., 42., 47., 55., 65., 75., 200.,
        15., 36., 38., 42., 47., 55., 65., 75., 200.,
        15., 36., 38., 42., 47., 55., 65., 75., 200.,
        15., 36., 38., 42., 47., 55., 65., 75., 200.,
        15., 36., 38., 42., 47., 55., 65., 75., 200.,
        15., 36., 38., 42., 47., 55., 65., 75., 200.,
        15., 36., 38., 42., 47., 55., 65., 75., 200.,
        15., 36., 38., 42., 47., 55., 65., 75., 200.,
        15., 36., 38., 42., 47., 55., 65., 75., 200.,
        15., 36., 38., 42., 47., 55., 65., 75., 200.,
        15., 36., 38., 42., 47., 55., 65., 75., 200.,
        ]

    print(TRIG.to_evaluator())
    print("-" * 120)
    print("nominal --->","eta:",etas,"r9:",rs,"pt:",pts,"result:",TRIG.to_evaluator().evaluate("nominal", etas, rs, pts))
    print("-" * 120)
    print("up      --->","eta:",etas,"r9:",rs,"pt:",pts,"result:",TRIG.to_evaluator().evaluate("up", etas, rs, pts))
    print("-" * 120)
    print("down    --->","eta:",etas,"r9:",rs,"pt:",pts,"result:",TRIG.to_evaluator().evaluate("down", etas, rs, pts))
    print("-" * 120)
    print()

    cset = cs.CorrectionSet(
        schema_version=2,
        description="Trigger SF",
        compound_corrections=None,
        corrections=[
            TRIG,
        ],
    )

    if opt.read == "":
        with open("TriggerSF_lead.json", "w") as fout:
            fout.write(cset.model_dump_json(exclude_unset=True))

        import gzip

        with gzip.open("TriggerSF_lead.json.gz", "wt") as fout:
            fout.write(cset.model_dump_json(exclude_unset=True))

    if opt.read != "":
        print("Reading back...")
        file = opt.read
        ceval = correctionlib.CorrectionSet.from_file(file)
        for corr in ceval.values():
            print(f"Correction {corr.name} has {len(corr.inputs)} inputs")
            for ix in corr.inputs:
                print(f"   Input {ix.name} ({ix.type}): {ix.description}")

        print("-" * 120)
        print("nominal --->","eta:",etas,"r9:",rs,"result:",ceval["TriggerSF"].evaluate("nominal", etas, rs, pts))
        print("-" * 120)
        print("up      --->","eta:",etas,"r9:",rs,"result:",ceval["TriggerSF"].evaluate("up", etas, rs, pts))
        print("-" * 120)
        print("down    --->","eta:",etas,"r9:",rs,"result:",ceval["TriggerSF"].evaluate("down", etas, rs, pts))
        print("-" * 120)


# trigger SF and uncertainty for UL2017 sublead photon. Dt:17/11/2020 fromflashgg: https://github.com/cms-analysis/flashgg/blob/4edea8897e2a4b0518dca76ba6c9909c20c40ae7/Systematics/python/flashggDiPhotonSystematics2017_Legacy_cfi.py#L27
# a lot of bins (81, 3 in SCeta, 3 in r9 and 9 in Pt)
# link to the presentation: https://indico.cern.ch/event/963617/contributions/4103623/attachments/2141570/3608645/Zee_Validation_UL2017_Update_09112020_Prasant.pdf

if opt.do_trigger_sublead or opt.do_trigger:
    inputs_ = ["SCEta", "r9", "pt"]
    edges_ = [[0.0, 1.5, 3, 999.0], [0.0, 0.56, 0.85, 0.9, 999.0], [0., 28., 31., 35., 40., 45., 50., 60., 70., 90., 999999.0]]
    content_ = {
        "nominal": [ #  35            37            40            45            50            60            70            90             inf    pt
            0.6155988939, 0.7165819087, 0.7381962831, 0.7671925006, 0.7999358222, 0.8254675016, 0.8297030540, 0.8451584417, 0.8522482004, 0.8871193652, # 0 < eta < 1.5, R9 < 0.56
            0.9028486970, 0.9739174387, 0.9756211698, 0.9785859435, 0.9814869681, 0.9836603606, 0.9808533747, 0.9788313651, 0.9766053770, 0.9667617117, # 0 < eta < 1.5, 0.56 < R9 < 0.85
            0.8933536854, 0.9973958724, 0.9980479262, 0.9987289490, 0.9992636424, 0.9994686970, 0.9995552559, 0.9992541003, 0.9996086647, 0.9996779894, # 0 < eta < 1.5, 0.85 < R9 < 0.90
            0.8933536854, 0.9973958724, 0.9980479262, 0.9987289490, 0.9992636424, 0.9994686970, 0.9995552559, 0.9992541003, 0.9996086647, 0.9996779894, # 0 < eta < 1.5, R9 > 0.90
            0.6100544113, 0.7427840769, 0.7761341323, 0.8117452882, 0.8319088440, 0.8583582498, 0.8736432627, 0.8907409748, 0.9046665266, 0.9190711276, # 1.5 < eta < 3, R9 < 0.56
            0.6100544113, 0.7427840769, 0.7761341323, 0.8117452882, 0.8319088440, 0.8583582498, 0.8736432627, 0.8907409748, 0.9046665266, 0.9190711276, # 1.5 < eta < 3, 0.56 < R9 < 0.85
            0.8283101205, 0.9538552575, 0.9597166341, 0.9617373097, 0.9624428298, 0.9581303007, 0.9621293579, 0.9670262230, 0.9721855102, 0.9753380476, # 1.5 < eta < 3, 0.85 < R9 < 0.90
            0.8582078363, 0.9911788518, 0.9961663139, 0.9974520554, 0.9983872590, 0.9988958563, 0.9987919975, 0.9992790060, 0.9994720350, 0.9995989436, # 1.5 < eta < 3, R9 > 0.90
            1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, # eta > 3, R9 < 0.56
            1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, # eta > 3, 0.56 < R9 < 0.85
            1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, # eta > 3, 0.85 < R9 < 0.90
            1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000  # eta > 3, R9 > 0.90
            ],
        "up": [ #  35            37            40            45            50            60            70            90             inf    pt
            0.6214049860, 0.7200945062, 0.7414391429, 0.7702228640, 0.8020985193, 0.8290374613, 0.8409052024, 0.8613078310, 0.8717261113, 0.9167756020, # 0 < eta < 1.5, R9 < 0.56
            0.9038765514, 0.9757861251, 0.9766339755, 0.9796238611, 0.9825145230, 0.9848254173, 0.9818546188, 0.9812502236, 0.9778944586, 0.9689814201, # 0 < eta < 1.5, 0.56 < R9 < 0.85
            0.8999503055, 0.9991625217, 0.9991339772, 0.9997289596, 1.0002637084, 1.0004687794, 1.0005552634, 1.0003072356, 1.0006087172, 1.0006779895, # 0 < eta < 1.5, 0.85 < R9 < 0.90
            0.8999503055, 0.9991625217, 0.9991339772, 0.9997289596, 1.0002637084, 1.0004687794, 1.0005552634, 1.0003072356, 1.0006087172, 1.0006779895, # 0 < eta < 1.5, R9 > 0.90
            0.6145174735, 0.7469034105, 0.7784538864, 0.8135725031, 0.8336459809, 0.8612610880, 0.8803866839, 0.8996338556, 0.9129986758, 0.9283278540, # 1.5 < eta < 3, R9 < 0.56
            0.6145174735, 0.7469034105, 0.7784538864, 0.8135725031, 0.8336459809, 0.8612610880, 0.8803866839, 0.8996338556, 0.9129986758, 0.9283278540, # 1.5 < eta < 3, 0.56 < R9 < 0.85
            0.8326968341, 0.9563037293, 0.9610060490, 0.9634283585, 0.9634432295, 0.9597145494, 0.9651548567, 0.9709326862, 0.9753379325, 0.9785131581, # 1.5 < eta < 3, 0.85 < R9 < 0.90
            0.8592510525, 0.9930106233, 0.9972930973, 0.9984525944, 0.9995239916, 0.9998959452, 0.9997922552, 1.0003196950, 1.0004723952, 1.0005989451, # 1.5 < eta < 3, R9 > 0.90
            1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, # eta > 3, R9 < 0.56
            1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, # eta > 3, 0.56 < R9 < 0.85
            1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, # eta > 3, 0.85 < R9 < 0.90
            1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000  # eta > 3, R9 > 0.90
        ],
        "down": [ #  35            37            40            45            50            60            70            90             inf    pt
            0.6097928018, 0.7130693112, 0.7349534233, 0.7641621372, 0.7977731251, 0.8218975419, 0.8185009056, 0.8290090524, 0.8327702895, 0.8574631284, # 0 < eta < 1.5, R9 < 0.56
            0.9018208426, 0.9720487523, 0.9746083641, 0.9775480259, 0.9804594132, 0.9824953039, 0.9798521306, 0.9764125066, 0.9753162954, 0.9645420033, # 0 < eta < 1.5, 0.56 < R9 < 0.85
            0.8867570653, 0.9956292231, 0.9969618752, 0.9977289384, 0.9982635764, 0.9984686146, 0.9985552484, 0.9982009650, 0.9986086122, 0.9986779893, # 0 < eta < 1.5, 0.85 < R9 < 0.90
            0.8867570653, 0.9956292231, 0.9969618752, 0.9977289384, 0.9982635764, 0.9984686146, 0.9985552484, 0.9982009650, 0.9986086122, 0.9986779893, # 0 < eta < 1.5, R9 > 0.90
            0.6055913491, 0.7386647433, 0.7738143782, 0.8099180733, 0.8301717071, 0.8554554116, 0.8668998415, 0.8818480940, 0.8963343774, 0.9098144012, # 1.5 < eta < 3, R9 < 0.56
            0.6055913491, 0.7386647433, 0.7738143782, 0.8099180733, 0.8301717071, 0.8554554116, 0.8668998415, 0.8818480940, 0.8963343774, 0.9098144012, # 1.5 < eta < 3, 0.56 < R9 < 0.85
            0.8239234069, 0.9514067857, 0.9584272192, 0.9600462609, 0.9614424301, 0.9565460520, 0.9591038591, 0.9631197598, 0.9690330879, 0.9721629371, # 1.5 < eta < 3, 0.85 < R9 < 0.90
            0.8571646201, 0.9893470803, 0.9950395305, 0.9964515164, 0.9972505264, 0.9978957674, 0.9977917398, 0.998238317, 0.9984716748, 0.9985989421, # 1.5 < eta < 3, R9 > 0.90
            1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, # eta > 3, R9 < 0.56
            1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, # eta > 3, 0.56 < R9 < 0.85
            1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, # eta > 3, 0.85 < R9 < 0.90
            1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, # eta > 3, R9 > 0.90
        ]
    }
    flow_ = "clamp"

    TRIG = cs.Correction(
        name="TriggerSF",
        description="Sublead Photon Trigger Scale Factor",
        generic_formulas=None,
        version=1,
        inputs=[
            cs.Variable(
                name="systematic", type="string", description="Systematic variation"
            ),
            cs.Variable(name="SCEta", type="real", description="Photon Super Cluster eta absolute value"),
            cs.Variable(
                name="r9",
                type="real",
                description="Photon full 5x5 R9, ratio E3x3/ERAW, where E3x3 is the energy sum of the 3 by 3 crystals surrounding the supercluster seed crystal and ERAW is the raw energy sum of the supercluster",
            ),
            cs.Variable(
                name="pt",
                type="real",
                description="Photon pt",
            ),
        ],
        output=cs.Variable(
            name="Wcorr",
            type="real",
            description="Multiplicative correction to event weight (per-photon)",
        ),
        data=cs.Category(
            nodetype="category",
            input="systematic",
            content=[
                cs.CategoryItem(
                    key="up",
                    value=multibinning(inputs_, edges_, content_["up"], flow_)
                ),
                cs.CategoryItem(
                    key="down",
                    value=multibinning(inputs_, edges_, content_["down"], flow_)
                ),
            ],
            default=multibinning(inputs_, edges_, content_["nominal"], flow_),
        ),
    )

    rich.print(TRIG)

    etas = [
        0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
        0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
        0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
        0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
        1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7,
        1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7,
        1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7,
        1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7,
        10., 10., 10., 10., 10., 10., 10., 10., 10., 10.,
        10., 10., 10., 10., 10., 10., 10., 10., 10., 10.,
        10., 10., 10., 10., 10., 10., 10., 10., 10., 10.,
        10., 10., 10., 10., 10., 10., 10., 10., 10., 10.
        ]
    rs = [
        0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
        0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
        0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88,
        0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99,
        0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
        0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
        0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88,
        0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99,
        0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
        0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
        0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88,
        0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99
        ]
    pts = [
        15., 30., 34., 38., 42., 47., 55., 65., 75., 200.,
        15., 30., 34., 38., 42., 47., 55., 65., 75., 200.,
        15., 30., 34., 38., 42., 47., 55., 65., 75., 200.,
        15., 30., 34., 38., 42., 47., 55., 65., 75., 200.,
        15., 30., 34., 38., 42., 47., 55., 65., 75., 200.,
        15., 30., 34., 38., 42., 47., 55., 65., 75., 200.,
        15., 30., 34., 38., 42., 47., 55., 65., 75., 200.,
        15., 30., 34., 38., 42., 47., 55., 65., 75., 200.,
        15., 30., 34., 38., 42., 47., 55., 65., 75., 200.,
        15., 30., 34., 38., 42., 47., 55., 65., 75., 200.,
        15., 30., 34., 38., 42., 47., 55., 65., 75., 200.,
        15., 30., 34., 38., 42., 47., 55., 65., 75., 200.,
        ]

    print(TRIG.to_evaluator())
    print("-" * 120)
    print("nominal --->","eta:",etas,"r9:",rs,"pt:",pts,"result:",TRIG.to_evaluator().evaluate("nominal", etas, rs, pts))
    print("-" * 120)
    print("up      --->","eta:",etas,"r9:",rs,"pt:",pts,"result:",TRIG.to_evaluator().evaluate("up", etas, rs, pts))
    print("-" * 120)
    print("down    --->","eta:",etas,"r9:",rs,"pt:",pts,"result:",TRIG.to_evaluator().evaluate("down", etas, rs, pts))
    print("-" * 120)
    print()

    cset = cs.CorrectionSet(
        schema_version=2,
        description="Trigger SF",
        compound_corrections=None,
        corrections=[
            TRIG,
        ],
    )

    if opt.read == "":
        with open("TriggerSF_sublead.json", "w") as fout:
            fout.write(cset.model_dump_json(exclude_unset=True))

        import gzip

        with gzip.open("TriggerSF_sublead.json.gz", "wt") as fout:
            fout.write(cset.model_dump_json(exclude_unset=True))

    if opt.read != "":
        print("Reading back...")
        file = opt.read
        ceval = correctionlib.CorrectionSet.from_file(file)
        for corr in ceval.values():
            print(f"Correction {corr.name} has {len(corr.inputs)} inputs")
            for ix in corr.inputs:
                print(f"   Input {ix.name} ({ix.type}): {ix.description}")

        print("-" * 120)
        print("nominal --->","eta:",etas,"r9:",rs,"result:",ceval["TriggerSF"].evaluate("nominal", etas, rs, pts))
        print("-" * 120)
        print("up      --->","eta:",etas,"r9:",rs,"result:",ceval["TriggerSF"].evaluate("up", etas, rs, pts))
        print("-" * 120)
        print("down    --->","eta:",etas,"r9:",rs,"result:",ceval["TriggerSF"].evaluate("down", etas, rs, pts))
        print("-" * 120)
