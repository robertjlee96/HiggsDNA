#!/usr/bin/env python
import argparse
import os
from higgs_dna.utils.logger_utils import setup_logger
import urllib.request
import pathlib
import shutil

parser = argparse.ArgumentParser(
    description="Simple utility script to retrieve the needed files for corections, luminostiy mask, systematics uncertainties ..."
)

parser.add_argument(
    "-t",
    "--target",
    dest="target",
    help="Choose the target to download (default: %(default)s)",
    default="GoldenJson",
    choices=["GoldenJSON", "cTag", "PhotonID", "PileupRun2", "SS"],
)

parser.add_argument(
    "-a",
    "--all",
    dest="all",
    action="store_true",
    help="Download all the targets (default: %(default)s)",
    default=False,
)
parser.add_argument(
    "--log", dest="log", type=str, default="INFO", help="Logger info level"
)
parser.add_argument(
    "--target_dir",
    type=str,
    default=None,
    help="directory to place the correction jsons, default: ../higgs-dna/systematics/JSONs",
)
parser.add_argument(
    "--analysis",
    type=str,
    default="higgs-dna-test",
    help="Name of the analysis you're perfoming, ideally it would match the output directory in which you're analysis parquet will end up, default: higgs-dna-test.",
)
parser.add_argument(
    "--log_dir",
    type=str,
    default="./json-log/",
    help="Log file summarising the json will end up here, default: ./json-log/",
)

args = parser.parse_args()

# ---------------------- A few helping functions  ----------------------


def fetch_file(target_name, logger, from_to_dict, type="url"):
    if type == "url":
        try:
            for ikey in from_to_dict.keys():
                with urllib.request.urlopen(from_to_dict[ikey]["from"]) as f:
                    json_object = f.read().decode("utf-8")
                # create the folder
                p = pathlib.Path(from_to_dict[ikey]["to"])
                p = pathlib.Path(*p.parts[:-1])  # remove file name
                p.mkdir(parents=True, exist_ok=True)
                with open(from_to_dict[ikey]["to"], "w") as f:
                    f.write(json_object)
                logger.info(
                    "[ {} ] {}: Download from {} to {}".format(
                        target_name,
                        ikey,
                        from_to_dict[ikey]["from"],
                        from_to_dict[ikey]["to"],
                    )
                )
        except:
            logger.info(
                "[ {} ] {}: Can't download from {}".format(
                    target_name, ikey, from_to_dict[ikey]["from"]
                )
            )
    elif type == "copy":
        try:
            for ikey in from_to_dict.keys():
                # create the folder
                p = pathlib.Path(from_to_dict[ikey]["to"])
                p = pathlib.Path(*p.parts[:-1])  # remove file name
                p.mkdir(parents=True, exist_ok=True)
                # copy
                shutil.copy(from_to_dict[ikey]["from"], from_to_dict[ikey]["to"])
                logger.info(
                    "[ {} ] {}: Copy from {} to {}".format(
                        target_name,
                        ikey,
                        from_to_dict[ikey]["from"],
                        from_to_dict[ikey]["to"],
                    )
                )
        except:
            logger.info(
                "[ {} ] {}: Can't copy from {}".format(
                    target_name, ikey, from_to_dict[ikey]["from"]
                )
            )


def get_ctag_json(logger, target_dir):
    if target_dir is not None:
        to_prefix = target_dir
    else:
        to_prefix = os.path.join(
            os.path.dirname(__file__), "../higgs_dna/systematics/JSONs"
        )

    from_to_dict = {
        "2016preVFP": {
            "from": "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2016preVFP_UL/ctagging.json.gz",
            "to": f"{to_prefix}/cTagSF/ctagging_2016preVFP.json.gz",
        },
        "2016postVFP": {
            "from": "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2016postVFP_UL/ctagging.json.gz",
            "to": f"{to_prefix}/cTagSF/ctagging_2016postVFP.json.gz",
        },
        "2017": {
            "from": "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2017_UL/ctagging.json.gz",
            "to": f"{to_prefix}/cTagSF/ctagging_2017.json.gz",
        },
        "2018": {
            "from": "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2018_UL/ctagging.json.gz",
            "to": f"{to_prefix}/cTagSF/ctagging_2018.json.gz",
        },
    }
    fetch_file("cTag", logger, from_to_dict, type="copy")


def get_photonid_json(logger, target_dir):
    if target_dir is not None:
        to_prefix = target_dir
    else:
        to_prefix = os.path.join(
            os.path.dirname(__file__), "../higgs_dna/systematics/JSONs"
        )

    from_to_dict = {
        "2017": {
            "from": "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2017_UL/photon.json.gz",
            "to": f"{to_prefix}/SF_photon_ID/photon.json.gz",
        },
    }
    fetch_file("PhotonID", logger, from_to_dict, type="copy")

def get_scale_and_smearing(logger, target_dir):
    if target_dir is not None:
        to_prefix = target_dir
    else:
        to_prefix = os.path.join(
            os.path.dirname(__file__), "../higgs_dna/systematics/scaleAndSmearing"
        )

    from_to_dict = {
        "2022FG": {
            "from": "/eos/cms/store/group/phys_egamma/akapoor/S+SJSON/Prompt2022FG/SS.json",
            "to": f"{to_prefix}/2022FG/SS.json",
        },
    }
    fetch_file("Scale and Smearing", logger, from_to_dict, type="copy")


def get_goldenjson(logger, target_dir):
    # References:
    # https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun3Analysis#Data
    # This is not really a correction JSON, so we only allow saving to a specific location
    # Commnenting out the code below, this was the previous method
    #if target_dir is not None:
    #    to_prefix = target_dir
    #else:
    #    to_prefix = os.path.join(
    #        os.path.dirname(__file__), "../metaconditions/pileup"
    #    )

    prefix = "../higgs_dna/metaconditions/CAF/certification/"
    from_to_dict = {
        "2022": {
            "from": "https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions22/Cert_Collisions2022_355100_362760_Golden.json",
            "to": os.path.join(
                prefix,
                "Collisions22/Cert_Collisions2022_355100_362760_Golden.json",
            ),
        },
        "2023": {
            "from": "https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions23/Cert_Collisions2023_366442_370092_Golden.json",
            "to": os.path.join(
                prefix,
                "Collisions23/Cert_Collisions2023_366442_370092_Golden.json",
            ),
        },
    }

    fetch_file("GoldenJSON", logger, from_to_dict, type="url")


def get_pileup_Run2(logger, target_dir):
    # Base URL for pileup JSONs
    base_path = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/LUM"

    if target_dir is not None:
        to_prefix = target_dir
    else:
        to_prefix = os.path.join(
            os.path.dirname(__file__), "../higgs_dna/metaconditions/pileup"
        )

    from_to_dict = {
        "2016preVFP": {
            "from": f"{base_path}/2016preVFP_UL/puWeights.json.gz",
            "to": f"{to_prefix}/pileup_2016preVFP.json.gz",
        },
        "2016postVFP": {
            "from": f"{base_path}/2016postVFP_UL/puWeights.json.gz",
            "to": f"{to_prefix}/pileup_2016postVFP.json.gz",
        },
        "2017": {
            "from": f"{base_path}/2017_UL/puWeights.json.gz",
            "to": f"{to_prefix}/pileup_2017.json.gz",
        },
        "2018": {
            "from": f"{base_path}/2018_UL/puWeights.json.gz",
            "to": f"{to_prefix}/pileup_2018.json.gz",
        },
    }

    fetch_file("Pileup_Run2", logger, from_to_dict, type="copy")


if __name__ == "__main__":
    # log output
    logfile = os.path.join(args.log_dir, f"{args.analysis}_jsons.log")
    p = pathlib.Path(logfile)
    p = pathlib.Path(*p.parts[:-1])  # remove file name
    p.mkdir(parents=True, exist_ok=True)

    logger = setup_logger(level=args.log, logfile=logfile, time=True)

    if args.all:
        get_goldenjson(logger, args.target_dir)
        get_pileup_Run2(logger, args.target_dir)
        get_scale_and_smearing(logger, args.target_dir)
        get_ctag_json(logger, args.target_dir)
        get_photonid_json(logger, args.target_dir)
        get_pileup_Run2(logger, args.target_dir)
    elif args.target == "GoldenJSON":
        get_goldenjson(logger, args.target_dir)
    elif args.target == "PileupRun2":
        get_pileup_Run2(logger, args.target_dir)
    elif args.target == "SS":
        get_scale_and_smearing(logger, args.target_dir)
    elif args.target == "cTag":
        get_ctag_json(logger, args.target_dir)
    elif args.target == "PhotonID":
        get_photonid_json(logger, args.target_dir)
    else:
        logger.info("Unknown target, exit now!")
        exit(0)

    logger.info(" " * 60)
    logger.info("Done")
    logger.info("-" * 60)


# example commands:
# python pull_files.py --all
# python pull_files.py --target GoldenJSON
# python pull_files.py --target cTag
# python pull_files.py --target GoldenJson --target_dir ./test_json --log_dir ./json-log --analysis goldenjson_test
