#!/usr/bin/env python
# Author Tiziano Bevilacqua (03/03/2023)
import os
import subprocess
from optparse import OptionParser
import json
from higgs_dna.utils.logger_utils import setup_logger
import inspect

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


# function to activate the (pre-existing) FlashggFinalFit environment and use Tree2WS scripts
# full Power Ranger style :D
def activate_final_fit(path, command):
    current_path = os.getcwd()
    os.chdir(path)
    os.system(
        f"eval `scram runtime -sh` && source {path}/flashggFinalFit/setup.sh && cd {path}/flashggFinalFit/Trees2WS && {command} "
    )
    os.chdir(current_path)


# --------------------------------------------------------------------------------------------------------------------------#
# - USAGE: -----------------------------------------------------------------------------------------------------------------#
# - python3 run_postprocess_steps.py --input ../out_dir_syst_090323/ --merged --root --ws --syst --cats --args "--do_syst" -#
# --------------------------------------------------------------------------------------------------------------------------#

# Read options from command line
usage = "Usage: python %prog filelists [options]"
parser = OptionParser(usage=usage)
parser.add_option("--input", dest="input", type="string", default="", help="input dir")
parser.add_option(
    "--merge",
    dest="merge",
    action="store_true",
    default=False,
    help="Do merging of the .parquet files",
)
parser.add_option(
    "--root",
    dest="root",
    action="store_true",
    default=False,
    help="Do root conversion step",
)
parser.add_option(
    "--ws",
    dest="ws",
    action="store_true",
    default=False,
    help="Do root to workspace conversion step",
)
parser.add_option(
    "--ws_config",
    dest="config",
    type="string",
    default="config_higgsdna_cats.py",
    help="configuration file for Tree2WS, as it is now it must be stored in Tree2WS directory in FinalFit",
)
parser.add_option(
    "--final_fit",
    dest="final_fit",
    type="string",
    default="/work/bevila_t/HpC_Analysis/test_flashggfinalfit/CMSSW_10_2_13/src",
    help="FlashggFinalFit path",
)  # the default is just for me, it should be changed but I don't see a way to make this generally valid
parser.add_option(
    "--syst",
    dest="syst",
    action="store_true",
    default=False,
    help="Do systematics variation treatment",
)
parser.add_option(
    "--cats",
    dest="cats",
    action="store_true",
    default=False,
    help="Split into categories",
)
parser.add_option(
    "--args",
    dest="args",
    type="string",
    default="",
    help="additional options for root converter: --do_syst, --notag",
)
parser.add_option(
    "--verbose",
    dest="verbose",
    type="string",
    default="INFO",
    help="verbose lefer for the logger: INFO (default), DEBUG",
)
(opt, args) = parser.parse_args()

if (opt.verbose != "INFO") and (opt.verbose != "DEBUG"):
    opt.verbose = "INFO"
logger = setup_logger(level=opt.verbose)

os.system(
    f"ls -l {opt.input} | tail -n +2 | grep -v .coffea | grep 201 | "
    + "awk '{print $NF}' > dirlist.txt"
)

EXEC_PATH = os.getcwd()
os.chdir(EXEC_PATH + "/" + opt.input)
IN_PATH = os.getcwd()
SCRIPT_DIR = os.path.dirname(
    os.path.abspath(inspect.getfile(inspect.currentframe()))
)  # script directory

process_dict = {
    "cH_4FS_FXFX_M125_2017": "ch",
    "ggh_M125_2017": "ggh",
    "tth_M125_2017": "tth",
    "vbf_M125_2017": "vbf",
    "vh_M125_2017": "vh",
}

var_dict = {
    "NOMINAL": "nominal",
    "FNUFUp": "FNUF_up",
    "FNUFDown": "FNUF_down",
    "ShowerShapeUp": "ShowerShape_up",
    "ShowerShapeDown": "ShowerShape_down",
}

# Here we prepare to split the output into categories, in the dictionary are defined the cuts to be applyed by pyarrow.ParquetDataset
# when reading the data, the variable obviously has to be in the dumped .parquet
# This can be improved by passing the configuration via json loading
if opt.cats:
    cat_dict = {
        "BDT_LOW_0J_PT_LOW": {
            "cat_filter": [
                ("bdt_score", ">", 0.8),
                ("bdt_score", "<=", 0.9),
                ("ge_one_jet_cut", "==", False),
                ("pt", "<=", 10),
            ]
        },
        "BDT_LOW_0J_PT_HIG": {
            "cat_filter": [
                ("bdt_score", ">", 0.8),
                ("bdt_score", "<=", 0.9),
                ("ge_one_jet_cut", "==", False),
                ("pt", ">", 10),
            ]
        },
        "BDT_HIG_0J_PT_LOW": {
            "cat_filter": [
                ("bdt_score", ">", 0.9),
                ("ge_one_jet_cut", "==", False),
                ("pt", "<=", 10),
            ]
        },
        "BDT_HIG_0J_PT_HIG": {
            "cat_filter": [
                ("bdt_score", ">", 0.9),
                ("ge_one_jet_cut", "==", False),
                ("pt", ">", 10),
            ]
        },
        "BDT_LOW_GE1J_PT_LOW": {
            "cat_filter": [
                ("bdt_score", ">", 0.8),
                ("bdt_score", "<=", 0.9),
                ("ge_one_jet_cut", "==", True),
                ("pt", "<=", 60),
            ]
        },
        "BDT_LOW_GE1J_PT_HIG": {
            "cat_filter": [
                ("bdt_score", ">", 0.8),
                ("bdt_score", "<=", 0.9),
                ("ge_one_jet_cut", "==", True),
                ("pt", ">", 60),
            ]
        },
        "BDT_HIG_GE1J_PT_LOW": {
            "cat_filter": [
                ("bdt_score", ">", 0.9),
                ("ge_one_jet_cut", "==", True),
                ("pt", "<=", 60),
            ]
        },
        "BDT_HIG_GE1J_PT_HIG": {
            "cat_filter": [
                ("bdt_score", ">", 0.9),
                ("ge_one_jet_cut", "==", True),
                ("pt", ">", 60),
            ]
        },
    }
else:
    cat_dict = {"NOTAG": {"cat_filter": [("pt", ">", -1.0)]}}
# I create a dictionary and save it to a temporary json so that this can be shared between the two scripts
# and then gets deleted to not leave trash around. We have to care for the environment :P.
# Not super elegant, open for suggestions
with open("category.json", "w") as file:
    file.write(json.dumps(cat_dict))

# Same for variation dictionary, which is shared between merging and ROOTing steps.
with open("variation.json", "w") as file:
    file.write(json.dumps(var_dict))

os.system(f"mv category.json {SCRIPT_DIR}/../higgs_dna/category.json")
os.system(f"mv variation.json {SCRIPT_DIR}/../higgs_dna/variation.json")
cat_dict = "category.json"

if opt.merge:
    with open(f"{EXEC_PATH}/dirlist.txt") as fl:
        files = fl.readlines()
        for file in files:
            file = file.split("\n")[0]
            if "M125" in file:
                if os.path.exists(f"{IN_PATH}/merged/{file}"):
                    raise Exception(
                        f"The selected target path: {IN_PATH}/merged/{file} already exists"
                    )

                MKDIRP(f"{IN_PATH}/merged/{file}")
                if opt.syst:
                    # if we have systematic variations in different files we have to split them in different directories
                    # otherwise they will be all merged at once in the same output file
                    for var in var_dict:
                        os.chdir(IN_PATH)
                        MKDIRP(f"{IN_PATH}/{file}/{var}")
                        os.system(
                            f"mv {IN_PATH}/{file}/*{var_dict[var]}.parquet {IN_PATH}/{file}/{var}"
                        )

                        MKDIRP(f"{IN_PATH}/merged/{file}/{var}")

                        os.chdir(SCRIPT_DIR)
                        os.system(
                            f"python3 merge_parquet.py --source {IN_PATH}/{file}/{var} --target {IN_PATH}/merged/{file}/{var}/ --cats {cat_dict}"
                        )

                else:
                    os.chdir(SCRIPT_DIR)
                    os.system(
                        f"python3 merge_parquet.py --source {IN_PATH}/{file} --target {IN_PATH}/merged/{file}/ --cats {cat_dict}"
                    )
            else:
                if os.path.exists(f"{IN_PATH}/merged/{file}/{file}_merged.parquet"):
                    raise Exception(
                        f"The selected target path: {IN_PATH}/merged/{file}/{file}_merged.parquet already exists"
                    )
                if not os.path.exists(f'{IN_PATH}/merged/Data_{file.split("_")[-1]}'):
                    MKDIRP(f'{IN_PATH}/merged/Data_{file.split("_")[-1]}')
                os.chdir(SCRIPT_DIR)
                os.system(
                    f'python3 merge_parquet.py --source {IN_PATH}/{file} --target {IN_PATH}/merged/Data_{file.split("_")[-1]}/{file}_ --cats {cat_dict}'
                )

        # at this point Data will be split in eras, here we merge them again in one allData file to rule them all
        os.system(
            f'python3 merge_parquet.py --source {IN_PATH}/merged/Data_{file.split("_")[-1]} --target {IN_PATH}/merged/Data_{file.split("_")[-1]}/allData_ --cats {cat_dict}'
        )

if opt.root:
    # Note, in my version of HiggsDNA I run the analysis splitting data per Era in different datasets
    # the treatment of data here is tested just with that structure
    with open(f"{EXEC_PATH}/dirlist.txt") as fl:
        files = fl.readlines()
        for file in files:
            file = file.split("\n")[0]
            if "M125" in file and file in process_dict:
                if os.path.exists(f"{IN_PATH}/root/{file}"):
                    raise Exception(
                        f"The selected target path: {IN_PATH}/root/{file} already exists"
                    )

                if os.listdir(f"{IN_PATH}/merged/{file}/"):
                    logger.info(f"Found merged files {IN_PATH}/merged/{file}/")
                else:
                    raise Exception(f"Merged parquet not found at {IN_PATH}/merged/")

                MKDIRP(f"{IN_PATH}/root/{file}")
                os.chdir(SCRIPT_DIR)
                os.system(
                    f"python3 convert_parquet_to_root.py {IN_PATH}/merged/{file}/merged.parquet {IN_PATH}/root/{file}/merged.root mc --process {process_dict[file]} {opt.args} --cats {cat_dict} --vars variation.json"
                )
            elif "M125" not in file:
                if os.listdir(f'{IN_PATH}/merged/Data_{file.split("_")[-1]}/'):
                    logger.info(
                        f'Found merged data files in: {IN_PATH}/merged/Data_{file.split("_")[-1]}/'
                    )
                else:
                    raise Exception(
                        f'Merged parquet not found at: {IN_PATH}/merged/Data_{file.split("_")[-1]}/'
                    )

                if os.path.exists(
                    f'{IN_PATH}/root/Data/allData_{file.split("_")[-1]}.root'
                ):
                    logger.info(
                        f'Data already converted: {IN_PATH}/root/Data/allData_{file.split("_")[-1]}.root'
                    )
                    continue
                elif not os.path.exists(f"{IN_PATH}/root/Data/"):
                    MKDIRP(f"{IN_PATH}/root/Data")
                    os.chdir(SCRIPT_DIR)
                    os.system(
                        f'python3 convert_parquet_to_root.py {IN_PATH}/merged/Data_{file.split("_")[-1]}/allData_merged.parquet {IN_PATH}/root/Data/allData_{file.split("_")[-1]}.root data --cats {cat_dict} --vars variation.json'
                    )
                else:
                    os.chdir(SCRIPT_DIR)
                    os.system(
                        f'python3 convert_parquet_to_root.py {IN_PATH}/merged/Data_{file.split("_")[-1]}/allData_merged.parquet {IN_PATH}/root/Data/allData_{file.split("_")[-1]}.root data --cats {cat_dict} --vars variation.json'
                    )

if opt.ws:
    if not os.listdir(opt.final_fit):
        raise Exception(
            f"The selected FlashggFinalFit path: {opt.final_fit} is invalid"
        )
    if os.path.exists(f"{IN_PATH}/root/Data"):
        os.system(f"echo Data >> {EXEC_PATH}/dirlist.txt")

    with open(f"{EXEC_PATH}/dirlist.txt") as fl:
        files = fl.readlines()
        if opt.syst:
            doSystematics = "--doSystematics"
        else:
            doSystematics = ""
        for dir in files:
            # skipping not important directories (they were used in previous steps)
            if "DoubleEG" in dir:
                logger.debug(f"Skipping directory {dir}")
                continue
            dir = dir.split("\n")[0]
            if "M125" in dir and dir in process_dict:
                if os.listdir(f"{IN_PATH}/root/{dir}/"):
                    filename = subprocess.check_output(
                        f"find {IN_PATH}/root/{dir} -name *.root -type f",
                        shell=True,
                        universal_newlines=True,
                    )
                else:
                    raise Exception(
                        f"The selected target path: {IN_PATH}/root/{dir} it's empty"
                    )
                command = f"python trees2ws.py --inputConfig {opt.config} --productionMode {process_dict[dir]} --year 2017 {doSystematics} --inputTreeFile {filename}"
                activate_final_fit(opt.final_fit, command)
            else:
                if os.listdir(f"{IN_PATH}/root/{dir}/"):
                    filename = subprocess.check_output(
                        f"find {IN_PATH}/root/{dir} -name *.root -type f",
                        shell=True,
                        universal_newlines=True,
                    )
                else:
                    raise Exception(
                        f"The selected target path: {IN_PATH}/root/{dir} it's empty"
                    )
                command = f"python trees2ws_data.py --inputConfig {opt.config} --inputTreeFile {filename}"
                activate_final_fit(opt.final_fit, command)
    os.chdir(EXEC_PATH)

# We don't want to leave trash around
if os.path.exists(f"{EXEC_PATH}/dirlist.txt"):
    os.system(f"rm {EXEC_PATH}/dirlist.txt")
if os.path.exists(f"{SCRIPT_DIR}/../higgs_dna/category.json"):
    os.system(f"rm {SCRIPT_DIR}/../higgs_dna/category.json")
if os.path.exists(f"{SCRIPT_DIR}/../higgs_dna/variation.json"):
    os.system(f"rm {SCRIPT_DIR}/../higgs_dna/variation.json")
