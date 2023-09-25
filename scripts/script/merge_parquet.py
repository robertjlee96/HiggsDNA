#!/usr/bin/env python
import argparse
import json
import os
import glob
import awkward
from higgs_dna.utils.logger_utils import setup_logger
import pyarrow.parquet as pq

parser = argparse.ArgumentParser(
    description="Simple utility script to merge all parquet files in one folder."
)
parser.add_argument(
    "--source",
    type=str,
    default="",
    help="Comma separated paths (with trailing slash) to folder where multiple parquet files are located. Careful: Folder should ONLY contain parquet files!",
)
parser.add_argument(
    "--target",
    type=str,
    default="",
    help="Comma separated paths (with trailing slash) to desired folder. Resulting merged file is placed there.",
)
parser.add_argument(
    "--cats",
    type=str,
    dest="cats_dict",
    default="",
    help="Dictionary containing category selections.",
)
parser.add_argument(
    "--is-data",
    default=False,
    action="store_true",
    help="Files to be merged are data and therefore do not require normalisation.",
)
parser.add_argument(
    "--skip-normalisation",
    default=False,
    action="store_true",
    help="Independent of file type, skip normalisation step",
)


args = parser.parse_args()
source_paths = args.source.split(",")
target_paths = args.target.split(",")

BASEDIR = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))
) + "/higgs_dna/"

logger = setup_logger(level="INFO")

if (
    (len(source_paths) != len(target_paths))
    or (args.source == "")
    or (args.target == "")
):
    logger.info("You gave a different number of sources and targets")
    exit


if args.cats_dict != "":
    with open(BASEDIR + "category.json") as pf:
        cat_dict = json.load(pf)
    for cat in cat_dict:
        logger.info(f"Found category: {cat}")
else:
    logger.info(
        "You provided an invalid dictionary containing categories information, have a look at your version of prepare_output_file.py"
    )
    logger.info(
        "An inclusive NOTAG category is used as default"
    )
    cat_dict = {"NOTAG": {"cat_filter": [("pt", ">", -1.0)]}}


# TODO: is it possible to read all files metadata with the ParquetDataset function. Currently extracting norm outside
if (not args.is_data) & (not args.skip_normalisation):
    logger.info(
        "Extracting sum of gen weights (before selection) from metadata of files to be merged."
    )
    sum_genw_beforesel_arr = []
    for i, source_path in enumerate(source_paths):
        source_files = glob.glob("%s/*.parquet" % source_path)
        sum_genw_beforesel = 0
        for f in source_files:
            sum_genw_beforesel += float(pq.read_table(f).schema.metadata[b'sum_genw_presel'])
        sum_genw_beforesel_arr.append(sum_genw_beforesel)
    logger.info(
        "Successfully extracted sum of gen weights (before selection)"
    )

for i, source_path in enumerate(source_paths):
    for cat in cat_dict:
        logger.info("-" * 125)
        logger.info(
            f"INFO: Starting parquet file merging. Attempting to read ParquetDataset from {source_path}, for category: {cat}"
        )
        dataset = pq.ParquetDataset(source_path, filters=cat_dict[cat]["cat_filter"])
        logger.info("ParquetDataset read successfully.")
        logger.info(
            f"Attempting to merge ParquetDataset and save to {target_paths[i]}."
        )
        pq.write_table(
            dataset.read(), target_paths[i] + cat + "_merged.parquet"
        )  # dataset.read() is a pyarrow table
        logger.info(
            f"Success! Merged parquet file is located in {target_paths[i]}{cat}_merged.parquet."
        )
        # If MC then open the merged dataset and add normalised weight column (sumw = efficiency)
        # TODO: can we add column before writing table and prevent re-reading in as awkward array
        if (not args.is_data) & (not args.skip_normalisation):
            # Remove ParquetDataset from memory and read file in as awkward array
            del dataset
            dataset_arr = awkward.from_parquet(target_paths[i] + cat + "_merged.parquet")
            # Add column for unnormalised weight
            dataset_arr['weight_nominal'] = dataset_arr['weight']
            dataset_arr['weight'] = dataset_arr['weight'] / sum_genw_beforesel_arr[i]
            awkward.to_parquet(dataset_arr, target_paths[i] + cat + "_merged.parquet")
            logger.info(
                "Successfully added normalised weight column to dataset"
            )
        logger.info("-" * 125)
