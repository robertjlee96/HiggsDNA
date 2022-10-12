import argparse
from higgs_dna.utils.logger_utils import setup_logger
import pandas as pd
import uproot

parser = argparse.ArgumentParser(description='Simple utility script to convert one parquet file into one ROOT file.')
parser.add_argument('source', type=str, help='Path to input file.')
parser.add_argument('target', type=str, help='Path to desired output file.')
args = parser.parse_args()
source_path = args.source
target_path = args.target

logger = setup_logger(level="INFO")

logger.info(f"Starting conversion of one parquet file to ROOT. Attempting to read file {source_path}.")
df = pd.read_parquet(source_path)
logger.info("Successfully read pandas dataframe from parquet file.")
dict = df.to_dict(orient='list') # to save it in proper format (default does it row-wise = bad)
logger.info("Successfully created dict from pandas dataframe.")
logger.info(f"Attempting to write dict to ROOT file {target_path}.")
with uproot.recreate(target_path) as file:
    file["myTree"] = dict
logger.info("Successfully converted parquet file to ROOT file.")
