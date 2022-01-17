from higgs_dna.utils.logger_utils import setup_logger

import argparse
import json
import os

xrootd_pfx = {
    "Americas": "root://cmsxrootd.fnal.gov/",
    "Eurasia": "root://xrootd-cms.infn.it/",
    "Yolo": "root://cms-xrd-global.cern.ch/",
}


def get_fetcher_args():
    parser = argparse.ArgumentParser(
        description="Query dasgoclient for dataset file lists"
    )

    parser.add_argument(
        "-i",
        "--input",
        help="What input dataset definition file to process.",
        required=True,
    )
    parser.add_argument(
        "-w",
        "--where",
        help="Where are you running your jobs? (default: %(default)s)",
        default="Americas",
        choices=["Americas", "Eurasia", "Yolo"],
    )
    parser.add_argument(
        "-x",
        "--xrootd",
        help="Override xrootd prefix with the one given.",
        default=None,
    )
    parser.add_argument(
        "--dbs-instance",
        dest="instance",
        help="The DBS instance to use for querying datasets. (default: %(default)s)",
        type=str,
        default="prod/global",
        choices=["prod/global", "prod/phys01", "prod/phys02", "prod/phys03"],
    )

    return parser.parse_args()


def get_dataset_dict(fset, xrd, dbs_instance):
    """
    Get a dictionary of dataset and the files in it.

    :param fset: A list of tuples with the format (dataset-short-name, path)
    :type fset: list
    :return fdict: A dictionary of dataset-short-name: list-of-files
    :rtype fdict: dict:
    """
    fdict = {}

    for name, dataset in fset:
        flist = (
            os.popen(
                (
                    "/cvmfs/cms.cern.ch/common/dasgoclient -query='instance={} file dataset={}'"
                ).format(dbs_instance, dataset)
            )
            .read()
            .split("\n")
        )
        if name not in fdict:
            fdict[name] = [xrd + f for f in flist if len(f) > 1]
        else:  # needed to collect all data samples into one common key "Data" (using append() would introduce a new element for the key)
            fdict[name].extend([xrd + f for f in flist if len(f) > 1])

    return fdict


if __name__ == "__main__":

    args = get_fetcher_args()

    logger = setup_logger(level="INFO")

    if ".txt" not in args.input:
        raise Exception("Input file must have '.txt' extension and be a text file!")

    fset = []
    with open(args.input) as fp:
        for i, line in enumerate(fp.readlines()):
            if line.strip().startswith("#"):
                continue
            fset.append(tuple(line.strip().split()))
            if len(fset[-1]) != 2:
                raise Exception(
                    f"Text file format should be '<short name> <dataset path>' and nothing else.\nInvalid spec on line {i+1}: '{line}'"
                )
    logger.info(f"Using following combination of datasets names and paths: {fset}")

    xrd = xrootd_pfx[args.where] if args.xrootd is None else args.xrootd
    logger.info(f"Using xrootd prefix: {xrd}")

    fdict = get_dataset_dict(fset, xrd, args.instance)

    # pprint.pprint(fdict, depth=1)
    with open(args.input[: args.input.rfind(".txt")] + ".json", "w") as fp:
        json.dump(fdict, fp, indent=4)
