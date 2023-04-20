import argparse
import json
from importlib import resources
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
args = parser.parse_args()
source_paths = args.source.split(",")
target_paths = args.target.split(",")

logger = setup_logger(level="INFO")

if (
    (len(source_paths) != len(target_paths))
    or (args.source == "")
    or (args.target == "")
):
    logger.info("You gave a different number of sources and targets")
    exit


if args.cats_dict != "":
    with resources.open_text("higgs_dna", "category.json") as pf:
        cat_dict = json.load(pf)
    for cat in cat_dict:
        logger.info(f"Found category: {cat}")
else:
    raise Exception(
        "You provided an invalid dictionary containing categories information, have a look at your version of prepare_output_file.py"
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
