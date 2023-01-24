from higgs_dna.workflows import workflows, taggers
from higgs_dna.metaconditions import metaconditions

import argparse
import sys
import os
import subprocess
import json
import re
import logging

logger = logging.getLogger(__name__)


def get_systematics_dict(string=None):
    """
    Function used by the runner to get which systematics affect which datasets.
    If string is a path to a json file, it is loaded.
    If string is a string, it must be in the form 'dataset1:systA,systB/dataset2:systC,systD'
    The dictionary returned is of the form {dataset1: [systA, systB], dataset2: [systC, systD]}

    :param string: string with systematics
    :type string: str
    :return: dictionary with systematics
    :rtype: dict
    """
    if string is None:
        return {}
    if string.endswith(".json"):
        with open(string) as f:
            return json.load(f)
    else:
        dct = {}
        dataset_column_systematics_list = re.split("/| ", string)
        for dataset_column_systematics in dataset_column_systematics_list:
            dataset_name, systematics_string = dataset_column_systematics.split(":")
            systematics_list = re.split(",|, ", systematics_string)
            dct[dataset_name] = systematics_list
        return dct


def get_main_parser():
    parser = argparse.ArgumentParser(
        description="Run Hgg Workflows on NanoAOD using coffea processor files"
    )
    # Analysis inputs
    parser.add_argument(
        "--json-analysis",
        dest="json_analysis_file",
        type=str,
        help="JSON analysis file where workflow, taggers, metaconditions, samples and systematics are defined.",
        default=None,
    )
    parser.add_argument(
        "--wf",
        "--workflow",
        dest="workflow",
        choices=list(workflows.keys()),
        help="Which processor to run",
        required=True if "--json-analysis" not in sys.argv else False,
    )
    parser.add_argument(
        "--ts",
        "--tagger-set",
        dest="taggers",
        nargs="+",
        default=None,
        choices=list(taggers.keys()),
        help="The tagger set to apply to this run.",
    )
    parser.add_argument(
        "--meta",
        "--metaconditions",
        dest="metaconditions",
        choices=list(metaconditions.keys()),
        help="What metaconditions to load",
        required=True if "--json-analysis" not in sys.argv else False,
    )
    parser.add_argument(
        "--samples",
        "--json",
        dest="samplejson",
        default="dummy_samples.json",
        help="JSON file containing dataset and file locations (default: %(default)s)",
    )
    parser.add_argument(
        "--systs",
        "--systematics",
        dest="systematics",
        default=None,
        type=str,
        help="Systematic variations file",
    )

    # File handling information
    parser.add_argument(
        "--no-trigger",
        dest="use_trigger",
        default=True,
        action="store_false",
        help="Turn off trigger selection",
    )
    parser.add_argument(
        "-d",
        "--dump",
        default=None,
        help="Path to dump parquet outputs to (default: None)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=r"output.coffea",
        help="Output filename (default: %(default)s)",
    )

    # Scale out
    parser.add_argument(
        "--executor",
        choices=[
            "iterative",
            "futures",
            "parsl/slurm",
            "parsl/condor",
            "dask/condor",
            "dask/slurm",
            "dask/lpc",
            "dask/lxplus",
            "dask/casa",  # Use for coffea-casa
            "vanilla_lxplus"
        ],
        default="futures",  # Local executor (named after concurrent futures package)
        help="The type of executor to use (default: %(default)s). Other options can be implemented. "
        "For example see https://parsl.readthedocs.io/en/stable/userguide/configuring.html"
        "- `parsl/slurm` - tested at DESY/Maxwell"
        "- `parsl/condor` - tested at DESY, RWTH"
        "- `dask/slurm` - tested at DESY/Maxwell"
        "- `dask/condor` - tested at DESY, RWTH"
        "- `dask/lpc` - custom lpc/condor setup (due to write access restrictions)"
        "- `dask/lxplus` - custom lxplus/condor setup (due to port restrictions)"
        "- `vanilla_lxplus` - custom plain lxplus submitter"
    )
    parser.add_argument(
        "-j",
        "--workers",
        type=int,
        default=12,
        help="Number of workers (cores/threads) to use for multi-worker executors "
        "(e.g. futures or condor) (default: %(default)s)",
    )
    parser.add_argument(
        "-m",
        "--memory",
        type=str,
        default="10GB",
        help="Memory to use for each job in distributed executors (default: %(default)s)",
    )
    parser.add_argument(
        "--walltime",
        type=str,
        default="01:00:00",
        help="Walltime to use for each job in distributed executors (default: %(default)s)",
    )
    parser.add_argument(
        "--disk",
        type=str,
        default="20GB",
        help="Disk space to use for each job in distributed executors (default: %(default)s)",
    )
    parser.add_argument(
        "-s",
        "--scaleout",
        type=int,
        default=6,
        help="Number of nodes to scale out to if using slurm/condor. Total number of "
        "concurrent threads is ``workers x scaleout`` (default: %(default)s)",
    )
    parser.add_argument(
        "--max-scaleout",
        dest="max_scaleout",
        type=int,
        default=250,
        help="The maximum number of nodes to adapt the cluster to. (default: %(default)s)",
    )
    parser.add_argument(
        "-q",
        "--queue",
        type=str,
        default=None,
        help="Queue to submit jobs to if using slurm/condor (default: %(default)s)",
    )
    parser.add_argument(
        "--voms",
        default=None,
        type=str,
        help="Path to voms proxy, accessible to worker nodes. Note that when this is specified "
        "the environment variable X509_CERT_DIR must be set to the certificates directory location",
    )

    # Debugging
    parser.add_argument(
        "--validate",
        action="store_true",
        default=False,
        help="Do not process, just check all files are accessible",
    )
    parser.add_argument("--skipbadfiles", action="store_true", help="Skip bad files.")
    parser.add_argument(
        "--only", type=str, default=None, help="Only process specific dataset or file"
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        metavar="N",
        help="Limit to the first N files of each dataset in sample JSON",
    )
    parser.add_argument(
        "--chunk",
        type=int,
        default=500000,
        metavar="N",
        help="Number of events per process chunk",
    )
    parser.add_argument(
        "--max",
        type=int,
        default=None,
        metavar="N",
        help="Max number of chunks to run in total",
    )
    parser.add_argument(
        "--skipCQR",
        default=False,
        action="store_true",
        help="Do not apply chained quantile regression (CQR) corrections",
    )
    parser.add_argument(
        "--debug",
        default=False,
        action="store_true",
        help="Print debug information with a logger",
    )

    return parser


def get_proxy():
    """
    Use voms-proxy-info to check if a proxy is available.
    If so, copy it to $HOME/.proxy and return the path.
    An exception is raised in the following cases:
    - voms-proxy-info is not installed
    - the proxy is not valid

    :return: Path to proxy
    :rtype: str
    """
    if subprocess.getstatusoutput("voms-proxy-info")[0] != 0:
        raise RuntimeError("voms-proxy-init not found. Please install it.")

    stat, out = subprocess.getstatusoutput("voms-proxy-info -e -p")
    # stat is 0 the proxy is valid
    if stat != 0:
        raise RuntimeError("No valid proxy found. Please create one.")

    _x509_localpath = out
    _x509_path = os.environ["HOME"] + f'/.{_x509_localpath.split("/")[-1]}'
    os.system(f"cp {_x509_localpath} {_x509_path}")

    logger.debug(f"Copied proxy to {_x509_path}")

    return _x509_path
