#!/usr/bin/env python
# Author Tiziano Bevilacqua (03/03/2023)
import os
import subprocess
from optparse import OptionParser
import json
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
# - python3 prepare_output_file.py --input ../out_dir_syst_090323/ --merge --root --ws --syst --cats --args "--do_syst" -#
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
    "--skip-normalisation",
    dest="skip_normalisation",
    action="store_true",
    default=False,
    help="Independent of file type, skip normalisation step",
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
    f"ls -l {opt.input} | tail -n +2 | grep -v .coffea | grep -v merged | grep -v root |"
    + "awk '{print $NF}' > dirlist.txt"
)

EXEC_PATH = os.getcwd()
os.chdir(opt.input)
IN_PATH = os.getcwd()
SCRIPT_DIR = os.path.dirname(
    os.path.abspath(__file__)
)  # script directory

process_dict = {
    "GluGluHtoGG_postEE_M125_2022": "ggh",
    "VBFHtoGG_postEE_M125_2022": "vbf",
    "ttHtoGG_postEE_M125_2022": "tth",
}

var_dict = {
    "NOMINAL": "nominal",
}

# Here we prepare to split the output into categories, in the dictionary are defined the cuts to be applyed by pyarrow.ParquetDataset
# when reading the data, the variable obviously has to be in the dumped .parquet
# This can be improved by passing the configuration via json loading
if opt.cats:
    cat_dict = {
        "EBEB_highR9highR9": {
            "cat_filter": [
                ("lead_isScEtaEB", "==", True),
                ("lead_r9", ">=", 0.85),
                ("sublead_isScEtaEB", "==", True),
                ("sublead_r9", ">=", 0.85),
            ]
        },
        "EBEB_highR9lowR9": {
            "cat_filter": [
                ("lead_isScEtaEB", "==", True),
                ("lead_r9", ">=", 0.85),
                ("sublead_isScEtaEB", "==", True),
                ("sublead_r9", "<", 0.85),
                ("sublead_r9", ">=", 0.5),
            ]
        },
        "EBEB_lowR9highR9": {
            "cat_filter": [
                ("lead_isScEtaEB", "==", True),
                ("lead_r9", "<", 0.85),
                ("lead_r9", ">=", 0.5),
                ("sublead_isScEtaEB", "==", True),
                ("sublead_r9", ">=", 0.85),
            ]
        },

        "EBEE_highR9highR9": {
            "cat_filter": [
                ("lead_isScEtaEB", "==", True),
                ("lead_r9", ">=", 0.85),
                ("sublead_isScEtaEE", "==", True),
                ("sublead_r9", ">=", 0.9),
            ]
        },
        "EBEE_highR9lowR9": {
            "cat_filter": [
                ("lead_isScEtaEB", "==", True),
                ("lead_r9", ">=", 0.85),
                ("sublead_isScEtaEE", "==", True),
                ("sublead_r9", "<", 0.9),
                ("sublead_r9", ">=", 0.8),
            ]
        },
        "EBEE_lowR9highR9": {
            "cat_filter": [
                ("lead_isScEtaEB", "==", True),
                ("lead_r9", "<", 0.85),
                ("lead_r9", ">=", 0.5),
                ("sublead_isScEtaEE", "==", True),
                ("sublead_r9", ">=", 0.9),
            ]
        },

        "EEEB_highR9highR9": {
            "cat_filter": [
                ("lead_isScEtaEE", "==", True),
                ("lead_r9", ">=", 0.9),
                ("sublead_isScEtaEB", "==", True),
                ("sublead_r9", ">=", 0.85),
            ]
        },
        "EEEB_highR9lowR9": {
            "cat_filter": [
                ("lead_isScEtaEE", "==", True),
                ("lead_r9", ">=", 0.9),
                ("sublead_isScEtaEB", "==", True),
                ("sublead_r9", "<", 0.85),
                ("sublead_r9", ">=", 0.5),
            ]
        },
        "EEEB_lowR9highR9": {
            "cat_filter": [
                ("lead_isScEtaEE", "==", True),
                ("lead_r9", "<", 0.9),
                ("lead_r9", ">=", 0.8),
                ("sublead_isScEtaEB", "==", True),
                ("sublead_r9", ">=", 0.85),
            ]
        },

        "EEEE_incl": {
            "cat_filter": [
                ("lead_isScEtaEE", "==", True),
                ("sublead_isScEtaEE", "==", True),
            ]
        }
    }
else:
    cat_dict = {"UNTAGGED": {"cat_filter": [("pt", ">", -1.0)]}}
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

# Define string if normalisation to be skipped
skip_normalisation_str = "--skip-normalisation" if opt.skip_normalisation else ""

if opt.merge:
    with open(f"{EXEC_PATH}/dirlist.txt") as fl:
        files = fl.readlines()
        for file in files:
            file = file.split("\n")[0]
            # MC dataset are supposed to have M125 in the name, this should be changed!
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

                        MKDIRP(f"{IN_PATH}/merged/{file}/{var_dict[var]}")

                        os.chdir(SCRIPT_DIR)
                        os.system(
                            f"python3 merge_parquet.py --source {IN_PATH}/{file}/{var_dict[var]} --target {IN_PATH}/merged/{file}/{var_dict[var]}/ --cats {cat_dict} {skip_normalisation_str}"
                        )

                else:
                    os.chdir(SCRIPT_DIR)
                    os.system(
                        f"python3 merge_parquet.py --source {IN_PATH}/{file}/nominal --target {IN_PATH}/merged/{file}/ --cats {cat_dict} {skip_normalisation_str}"
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
                    f'python3 merge_parquet.py --source {IN_PATH}/{file}/nominal --target {IN_PATH}/merged/Data_{file.split("_")[-1]}/{file}_ --cats {cat_dict} --is-data'
                )

        # at this point Data will be split in eras if any Data dataset is present, here we merge them again in one allData file to rule them all
        # we also skip this step if there is no Data
        for file in files:
            file = file.split("\n")[0]  # otherwise it contains an end of line and messes up the os.walk() call
            if "Data" in file or "DoubleEG" in file:
                dirpath, dirnames, filenames = next(os.walk(f'{IN_PATH}/merged/Data_{file.split("_")[-1]}'))
                if len(filenames) > 0:
                    os.system(
                        f'python3 merge_parquet.py --source {IN_PATH}/merged/Data_{file.split("_")[-1]} --target {IN_PATH}/merged/Data_{file.split("_")[-1]}/allData_ --cats {cat_dict} --is-data'
                    )
                    break
                else:
                    logger.info(f'No merged parquet found for {file} in the directory: {IN_PATH}/merged/Data_{file.split("_")[-1]}')

if opt.root:
    logger.info("Starting root step")
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

    data_done = False

    with open(f"{EXEC_PATH}/dirlist.txt") as fl:
        files = fl.readlines()
        if opt.syst:
            doSystematics = "--doSystematics"
        else:
            doSystematics = ""
        for dir in files:
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
            elif "Data" in dir and not data_done:
                if os.listdir(f"{IN_PATH}/root/Data/"):
                    filename = subprocess.check_output(
                        f"find {IN_PATH}/root/Data -name *.root -type f",
                        shell=True,
                        universal_newlines=True,
                    )
                else:
                    raise Exception(
                        f"The selected target path: {IN_PATH}/root/{dir} it's empty"
                    )
                command = f"python trees2ws_data.py --inputConfig {opt.config} --inputTreeFile {filename}"
                activate_final_fit(opt.final_fit, command)
                data_done = True
    os.chdir(EXEC_PATH)

# We don't want to leave trash around
if os.path.exists(f"{EXEC_PATH}/dirlist.txt"):
    os.system(f"rm {EXEC_PATH}/dirlist.txt")
if os.path.exists(f"{SCRIPT_DIR}/../higgs_dna/category.json"):
    os.system(f"rm {SCRIPT_DIR}/../higgs_dna/category.json")
if os.path.exists(f"{SCRIPT_DIR}/../higgs_dna/variation.json"):
    os.system(f"rm {SCRIPT_DIR}/../higgs_dna/variation.json")
