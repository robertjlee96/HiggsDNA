#!/usr/bin/env python
import argparse
import os
import sys
from higgs_dna.utils.logger_utils import setup_logger
import urllib.request
import pathlib
import shutil
import subprocess
from distutils.dir_util import copy_tree


parser = argparse.ArgumentParser(
    description="Simple utility script to retrieve the needed files for corections, luminostiy mask, systematics uncertainties ..."
)

parser.add_argument(
    "-t",
    "--target",
    dest="target",
    help="Choose the target to download (default: %(default)s)",
    default="GoldenJson",
    choices=["GoldenJSON", "cTag", "PhotonID", "PU_run2", "SS", "JetMET","CDFs", "JEC", "JER", "Material", "TriggerSF", "PreselSF", "eVetoSF", "Flows"],
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
    help="Log file summarising the json will end up here, default: ./json-log/"
)

args = parser.parse_args()


# ---------------------- A few helping functions  ----------------------
def unzip_gz_with_zcat(logger, input_file, output_file):
    try:
        # Check if zcat is available in the system
        subprocess.check_call(["zcat", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # Run zcat command to unzip the gz file
        with open(output_file, "wb") as output:
            subprocess.check_call(["zcat", input_file], stdout=output)
        logger.info(f"File '{input_file}' successfully unzipped to '{output_file}'.")
        # Remove the gz file after extraction
        os.remove(input_file)
        logger.info(f"File '{input_file}' deleted.")
    except subprocess.CalledProcessError as e:
        logger.info(f"Error: {e}")
        sys.exit(1)


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
                if os.path.isdir(from_to_dict[ikey]["from"]):
                    copy_tree(from_to_dict[ikey]["from"], from_to_dict[ikey]["to"])
                else:
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


def get_jec_files(logger, target_dir):
    if target_dir is not None:
        to_prefix = target_dir
    else:
        to_prefix = os.path.join(
            os.path.dirname(__file__), "../higgs_dna/systematics/data/"
        )

    from_to_dict = {
        "2017": {
            "from": "/eos/cms/store/group/phys_higgs/cmshgg/tbevilac/JECDatabase/textFiles/Summer19UL17_V5_MC",
            "to": f"{to_prefix}/Summer19UL17_MC/JEC/",
        },
        "2022postEE": {
            "from": "/eos/cms/store/group/phys_higgs/cmshgg/tbevilac/JECDatabase/textFiles/Winter22Run3_V2_MC",
            "to": f"{to_prefix}/Winter22Run3_MC/JEC/",
        },
    }
    fetch_file("JEC", logger, from_to_dict, type="copy")
    os.system(f"rename .txt .junc.txt {to_prefix}/*/JEC/*/*Uncertainty*.txt")


def get_jer_files(logger, target_dir):
    if target_dir is not None:
        to_prefix = target_dir
    else:
        to_prefix = os.path.join(
            os.path.dirname(__file__), "../higgs_dna/systematics/data/"
        )

    from_to_dict = {
        "2017": {
            "from": "/eos/cms/store/group/phys_higgs/cmshgg/tbevilac/JRDatabase/textFiles/Summer19UL17_JRV2_MC",
            "to": f"{to_prefix}/Summer19UL17_MC/JER/",
        },
        "2022postEE": {
            "from": "/eos/cms/store/group/phys_higgs/cmshgg/tbevilac/JRDatabase/textFiles/JR_Winter22Run3_V1_MC",
            "to": f"{to_prefix}/Winter22Run3_MC/JER/",
        },
    }
    fetch_file("JER", logger, from_to_dict, type="copy")
    os.system(f"rename PtResolution_AK4PFchs.txt PtResolution_AK4PFchs.jr.txt {to_prefix}/*/JER/*/*PtResolution_AK4PFchs.txt")
    os.system(f"rename PtResolution_AK4PFPuppi.txt PtResolution_AK4PFPuppi.jr.txt {to_prefix}/*/JER/*/*PtResolution_AK4PFPuppi.txt")
    os.system(f"rename SF_AK4PFchs.txt SF_AK4PFchs.jersf.txt {to_prefix}/*/JER/*/*SF_AK4PFchs.txt")
    os.system(f"rename SF_AK4PFPuppi.txt SF_AK4PFPuppi.jersf.txt {to_prefix}/*/JER/*/*SF_AK4PFPuppi.txt")


def get_material_json(logger, target_dir):
    if target_dir is not None:
        to_prefix = target_dir
    else:
        to_prefix = os.path.join(
            os.path.dirname(__file__), "../higgs_dna/systematics/JSONs/"
        )

    from_to_dict = {
        "2017": {
            "from": "/eos/cms/store/group/phys_higgs/cmshgg/tbevilac/JSONs/Material.json",
            "to": f"{to_prefix}",
        }
    }
    fetch_file("Material", logger, from_to_dict, type="copy")


def get_trigger_json(logger, target_dir):
    if target_dir is not None:
        to_prefix = target_dir
    else:
        to_prefix = os.path.join(
            os.path.dirname(__file__), "../higgs_dna/systematics/JSONs/"
        )

    from_to_dict = {
        "2017": {
            "from": "/eos/cms/store/group/phys_higgs/cmshgg/tbevilac/JSONs/TriggerSF*.json",
            "to": f"{to_prefix}",
        }
    }
    fetch_file("TriggerSF", logger, from_to_dict, type="copy")


def get_presel_json(logger, target_dir):
    if target_dir is not None:
        to_prefix = target_dir
    else:
        to_prefix = os.path.join(
            os.path.dirname(__file__), "../higgs_dna/systematics/JSONs/"
        )

    from_to_dict = {
        "2017": {
            "from": "/eos/cms/store/group/phys_higgs/cmshgg/tbevilac/JSONs/PreselSF.json",
            "to": f"{to_prefix}",
        }
    }
    fetch_file("PreselSF", logger, from_to_dict, type="copy")


def get_eveto_json(logger, target_dir):
    if target_dir is not None:
        to_prefix = target_dir
    else:
        to_prefix = os.path.join(
            os.path.dirname(__file__), "../higgs_dna/systematics/JSONs/"
        )

    from_to_dict = {
        "2017": {
            "from": "/eos/cms/store/group/phys_higgs/cmshgg/tbevilac/JSONs/eVetoSF.json",
            "to": f"{to_prefix}",
        }
    }
    fetch_file("eVetoSF", logger, from_to_dict, type="copy")


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
    # see https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammSFandSSRun3#Scale_And_Smearings_Correctionli
    if target_dir is not None:
        to_prefix = target_dir
    else:
        to_prefix = os.path.join(
            os.path.dirname(__file__), "../higgs_dna/systematics/JSONs/scaleAndSmearing"
        )

    from_to_dict = {
        "2022preEE": {
            "from": "/eos/cms/store/group/phys_egamma/akapoor/S+SJSON/2022Re-recoBCD/photonSS.json.gz",
            "to": f"{to_prefix}/SS_Rereco2022BCD.json.gz",
        },
        "2022postEE": {
            "from": "/eos/cms/store/group/phys_egamma/akapoor/S+SJSON/2022Re-recoE+PromptFG/photonSS.json.gz",
            "to": f"{to_prefix}/SS_RerecoE_PromptFG_2022.json.gz",
        },
    }
    fetch_file("Scale and Smearing", logger, from_to_dict, type="copy")
    # Now, unpack the gz to have the raw JSONs
    unzip_gz_with_zcat(logger, f"{to_prefix}/SS_Rereco2022BCD.json.gz", f"{to_prefix}/SS_Rereco2022BCD.json")
    unzip_gz_with_zcat(logger, f"{to_prefix}/SS_RerecoE_PromptFG_2022.json.gz", f"{to_prefix}/SS_RerecoE_PromptFG_2022.json")


def get_mass_decorrelation_CDF(logger, target_dir):
    if target_dir is not None:
        to_prefix = target_dir
    else:
        to_prefix = os.path.join(
            os.path.dirname(__file__), "../higgs_dna/tools"
        )

    from_to_dict = {
        "2022FG": {
            "from": "/eos/cms/store/group/phys_higgs/cmshgg/earlyRun3Hgg/mass_decorrelation/CDFs/Smeared_Diphoton_CDFs.pkl.gz",
            "to": f"{to_prefix}/",
        },
    }
    fetch_file("CDFs", logger, from_to_dict, type="copy")


def get_Flow_files(logger, target_dir):
    if target_dir is not None:
        to_prefix = target_dir
    else:
        to_prefix = os.path.join(
            os.path.dirname(__file__), "../higgs_dna/tools/flows"
        )

    from_to_dict = {
        "2022FG": {
            "from": "/eos/cms/store/group/phys_higgs/cmshgg/earlyRun3Hgg/flow_corrections/",
            "to": f"{to_prefix}/",
        },
    }
    fetch_file("Flows", logger, from_to_dict, type="copy")

def get_goldenjson(logger, target_dir):
    # References:
    # https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun3Analysis#Data
    # This is not really a correction JSON, so we only allow saving to a specific location
    # Commnenting out the code below, this was the previous method
    # if target_dir is not None:
    #    to_prefix = target_dir
    # else:
    #    to_prefix = os.path.join(
    #        os.path.dirname(__file__), "../metaconditions/pileup"
    #    )

    prefix = "../higgs_dna/metaconditions/CAF/certification/"

    from_to_dict = {
        "2016": {
            "from": "https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt",
            "to": os.path.join(
                prefix,
                "Collisions16/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt",
            ),
        },
        "2017": {
            "from": "https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt",
            "to": os.path.join(
                prefix,
                "Collisions17/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt",
            ),
        },
        "2018": {
            "from": "https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt",
            "to": os.path.join(
                prefix,
                "Collisions18/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt",
            ),
        },
        "2022": {
            "from": "https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions22/Cert_Collisions2022_355100_362760_Golden.json",
            "to": os.path.join(
                prefix,
                "Collisions22/Cert_Collisions2022_355100_362760_Golden.json",
            ),
        },
        "2023": {
            "from": "https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions23/Cert_Collisions2023_366442_370790_Golden.json",
            "to": os.path.join(
                prefix,
                "Collisions23/Cert_Collisions2023_366442_370790_Golden.json",
            ),
        },
    }

    fetch_file("GoldenJSON", logger, from_to_dict, type="url")


def get_jetmet_json(logger, target_dir):
    # References:
    # json pog of JME: https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/JME
    # jetmapveto: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis#From_JME
    base_path = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME"
    if target_dir is not None:
        to_prefix = target_dir
    else:
        to_prefix = os.path.dirname(__file__)

    from_to_dict = {
        "2016preVFP": {
            "from": os.path.join(base_path, "2016preVFP_UL"),
            "to": os.path.join(
                to_prefix,
                "../higgs_dna/systematics/JSONs/POG/JME/2016preVFP_UL",
            ),
        },
        "2016postVFP": {
            "from": os.path.join(base_path, "2016postVFP_UL"),
            "to": os.path.join(
                to_prefix,
                "../higgs_dna/systematics/JSONs/POG/JME/2016postVFP_UL",
            ),
        },
        "2017": {
            "from": os.path.join(base_path, "2017_UL"),
            "to": os.path.join(
                to_prefix,
                "../higgs_dna/systematics/JSONs/POG/JME/2017_UL",
            ),
        },
        "2018": {
            "from": os.path.join(base_path, "2018_UL"),
            "to": os.path.join(
                to_prefix,
                "../higgs_dna/systematics/JSONs/POG/JME/2018_UL",
            ),
        },
        "2022Prompt": {
            "from": os.path.join(base_path, "2022_Prompt"),
            "to": os.path.join(
                to_prefix,
                "../higgs_dna/systematics/JSONs/POG/JME/2022_Prompt",
            ),
        },
    }

    fetch_file("JetMET", logger, from_to_dict, type="copy")


def get_pileup_Run2(logger, target_dir):
    # Base URL for pileup JSONs
    base_path = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/LUM"

    if target_dir is not None:
        to_prefix = target_dir
    else:
        to_prefix = os.path.join(
            os.path.dirname(__file__), "../higgs_dna/systematics/JSONs/pileup/"
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
        get_mass_decorrelation_CDF(logger, args.target_dir)
        get_Flow_files(logger, args.target_dir)
        get_ctag_json(logger, args.target_dir)
        get_photonid_json(logger, args.target_dir)
        get_pileup_Run2(logger, args.target_dir)
        get_jetmet_json(logger, args.target_dir)
    elif args.target == "GoldenJSON":
        get_goldenjson(logger, args.target_dir)
    elif args.target == "PU_run2":
        get_pileup_Run2(logger, args.target_dir)
    elif args.target == "SS":
        get_scale_and_smearing(logger, args.target_dir)
    elif args.target == "CDFs":
        get_mass_decorrelation_CDF(logger, args.target_dir)
    elif args.target == "Flows":
        get_Flow_files(logger, args.target_dir)
    elif args.target == "cTag":
        get_ctag_json(logger, args.target_dir)
    elif args.target == "PhotonID":
        get_photonid_json(logger, args.target_dir)
    elif args.target == "JetMET":
        get_jetmet_json(logger, args.target_dir)
    elif args.target == "JEC":
        get_jec_files(logger, args.target_dir)
    elif args.target == "JER":
        get_jer_files(logger, args.target_dir)
    elif args.target == "Material":
        get_material_json(logger, args.target_dir)
    elif args.target == "TriggerSF":
        get_trigger_json(logger, args.target_dir)
    elif args.target == "PreselSF":
        get_presel_json(logger, args.target_dir)
    elif args.target == "eVetoSF":
        get_eveto_json(logger, args.target_dir)
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
