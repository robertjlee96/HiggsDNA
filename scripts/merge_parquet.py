import argparse
from higgs_dna.utils.logger_utils import setup_logger
import pyarrow.parquet as pq

parser = argparse.ArgumentParser(description='Simple utility script to merge all parquet files in one folder.')
parser.add_argument('source', type=str, help='Path (with trailing slash) to folder where multiple parquet files are located. Careful: Folder should ONLY contain parquet files!')
parser.add_argument('target', type=str, help='Path (with trailing slash) to desired folder. Resulting merged file is placed there.')
args = parser.parse_args()
source_path = args.source
target_path = args.target

logger = setup_logger(level="INFO")

logger.info(f"INFO: Starting parquet file merging. Attempting to read ParquetDataset from {source_path}.")
dataset = pq.ParquetDataset(source_path)
logger.info("ParquetDataset read successfully.")
logger.info(f"Attempting to merge ParquetDataset and save to {target_path}.")
pq.write_table(dataset.read(), target_path+'merged.parquet') # dataset.read() is a pyarrow table
logger.info(f"Success! Merged parquet file is located in {target_path}.")