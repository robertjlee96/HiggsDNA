# import root_pandas
import numpy as np
import pandas as pd
from decorrelator import cdfCalc
import argparse
import glob


def main(options):

    df = pd.DataFrame()

    files = glob.glob(options.infile)
    data = [pd.read_parquet(f) for f in files]
    events = pd.concat(data,ignore_index=True)

    df["sigma_m_over_m"] = events[options].to_numpy()
    df["mass"] = events.mass.to_numpy()
    df["weight"] = events.weight.to_numpy()

    print(f"INFO: found {len(events)} events")

    df["sigma_m_over_m"] = events[options].to_numpy()

    df["mass"] = events.mass.to_numpy()
    df["weight"] = events.weight.to_numpy()

    # Evaluating and dumping the CDFs in bins of mass
    calc = cdfCalc(df, options.tree,'mass',np.linspace(100, 180, 161))
    calc.dumpCdfs(options.cdfsFile)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    requiredArgs = parser.add_argument_group()
    requiredArgs.add_argument('-i', '--infile', action='store', type=str, required=True)
    requiredArgs.add_argument('-i', '--infile', action='store', type=str, required=True)
    requiredArgs.add_argument('-c', '--cdfsFile', action='store', type=str, required=True)
    requiredArgs.add_argument('-t', '--tree', action='store', type=str, required=True)
    options = parser.parse_args()
    main(options)
