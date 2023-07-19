#!/usr/bin/env python
import argparse
import datetime
import os
from higgs_dna.utils.logger_utils import setup_logger

# ---------------------- A few helping functions  ----------------------

def MKDIRP(dirpath, verbose=False, dry_run=False):
    if verbose:
        print("\033[1m" + ">" + "\033[0m" + ' os.mkdirs("' + dirpath + '")')
    if dry_run:
        return
    try:
        os.makedirs(dirpath)
    except OSError:
        if not os.path.isdir(dirpath):
            raise
    return

parser = argparse.ArgumentParser(
    description="Simple utility script to retrieve the needed json files for the systematics uncertainties."
)
parser.add_argument(
    "--target_dir",
    type=str,
    default="../higgs_dna/systematics/JSONs",
    help="directory to place the correction jsons, default: ../higgs-dna/systematics/JSONs",
)
parser.add_argument(
    "--analysis",
    type=str,
    default="higgs-dna-test",
    help="Name of the analysis you're perfoming, ideally it would match the output directory in which you're analysis parquet will end up, default: higgs-dna-test.",
)
parser.add_argument(
    "--log",
    type=str,
    default="./json-log/",
    help="Log file summarising the json will end up here, default: ./json-log/",
)

args = parser.parse_args()

logger = setup_logger(level="INFO")
logger.info("-" * 60)
logger.info("{:^60}".format("Fetching jsons"))
logger.info("   _____                 . . . . . o o o o o")
logger.info("  __|[_]|__ ___________ _______    ____      o")
logger.info(" |[] [] []| [] [] [] [] [_____(__  ][]]_n_n__][.")
logger.info("_|________|_[_________]_[________]_|__|________)<")
logger.info("  oo    oo 'oo      oo ' oo    oo 'oo 0000---oo\_")
logger.info("~"*60)

logger.info("-" * 60)

json_dict = {
    "cTagSF": {
        "2016preVFP": {
            "from": "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2016preVFP_UL/ctagging.json.gz",
            "to": f"{args.target_dir}/cTagSF/ctagging_2016preVFP.json.gz"
        },
        "2016postVFP": {
            "from": "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2016postVFP_UL/ctagging.json.gz",
            "to": f"{args.target_dir}/cTagSF/ctagging_2016postVFP.json.gz"
        },
        "2017": {
            "from": "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2017_UL/ctagging.json.gz",
            "to": f"{args.target_dir}/cTagSF/ctagging_2017.json.gz"
        },
        "2018": {
            "from": "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2018_UL/ctagging.json.gz",
            "to": f"{args.target_dir}/cTagSF/ctagging_2018.json.gz"
        }
    },
    "SF_photon_ID": {
        "2017": {
            "from": "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2017_UL/photon.json.gz",
            "to": f"{args.target_dir}/SF_photon_ID/photon.json.gz"
        }
    },
}

if not os.path.exists(args.log):
    MKDIRP(args.log)
    logger.info(f"created log dir at: {args.log}")
else:
    logger.info(f"logs will be written here: {args.log}")

os.system(f"echo '|-----------------------------------------------------------------------|' > {args.log}/{args.analysis}_jsons.log")
os.system(f"echo '      Summary of json conditions for analysis: {args.analysis} ' >> {args.log}/{args.analysis}_jsons.log")
os.system(f"echo '|-----------------------------------------------------------------------|' >> {args.log}/{args.analysis}_jsons.log")
os.system(f"echo '' >> {args.log}/{args.analysis}_jsons.log")
os.system(f"echo 'this file was created on: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}' >> {args.log}/{args.analysis}_jsons.log")
os.system(f"echo '' >> {args.log}/{args.analysis}_jsons.log")

for syst in json_dict:
    logger.info(f"looking for {syst} jsons...")
    os.system(f"echo 'looking for {syst} jsons...' >> {args.log}/{args.analysis}_jsons.log")
    for era in json_dict[syst]:
        if os.path.exists(json_dict[syst][era]["to"]):
            logger.info(f'The selected target path: {json_dict[syst][era]["to"]} already exists, skipping...')
            os.system(f"echo 'The selected target path: {json_dict[syst][era]['to']} already exists, skipping...' >> {args.log}/{args.analysis}_jsons.log")
            continue
        elif not os.path.exists(os.path.dirname(os.path.abspath(json_dict[syst][era]["to"]))):
            MKDIRP(os.path.dirname(os.path.abspath(json_dict[syst][era]["to"])))

        os.system(f'scp {json_dict[syst][era]["from"]} {json_dict[syst][era]["to"]}')

        logger.info(f'* found {json_dict[syst][era]["from"]} for era: {era}')
        logger.info(f'  ===> placed in: {json_dict[syst][era]["to"]}')
        os.system(f"echo '* found {json_dict[syst][era]['from']} for era: {era}' >> {args.log}/{args.analysis}_jsons.log")
        os.system(f"echo '  ===> placed in: {json_dict[syst][era]['to']}' >> {args.log}/{args.analysis}_jsons.log")

logger.info(" " * 60)
logger.info("Done") 
logger.info("-" * 60)