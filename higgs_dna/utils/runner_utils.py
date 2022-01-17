from higgs_dna.workflows import workflows, taggers
from higgs_dna.metaconditions import metaconditions
import argparse


def get_main_parser():
    parser = argparse.ArgumentParser(
        description="Run Hgg Workflows on NanoAOD using coffea processor files"
    )
    # Inputs
    parser.add_argument(
        "--wf",
        "--workflow",
        dest="workflow",
        choices=list(workflows.keys()),
        help="Which processor to run",
        required=True,
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
        required=True,
    )
    parser.add_argument(
        "--systs",
        "--systematics",
        dest="systematics",
        default=None,
        type=str,
        help="Systematic variations file",
    )
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
    parser.add_argument(
        "--samples",
        "--json",
        dest="samplejson",
        default="dummy_samples.json",
        help="JSON file containing dataset and file locations (default: %(default)s)",
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
        ],
        default="futures",  # Local executor (named after concurrent futures package)
        help="The type of executor to use (default: %(default)s). Other options can be implemented. "
        "For example see https://parsl.readthedocs.io/en/stable/userguide/configuring.html"
        "- `parsl/slurm` - tested at DESY/Maxwell"
        "- `parsl/condor` - tested at DESY, RWTH"
        "- `dask/slurm` - tested at DESY/Maxwell"
        "- `dask/condor` - tested at DESY, RWTH"
        "- `dask/lpc` - custom lpc/condor setup (due to write access restrictions)"
        "- `dask/lxplus` - custom lxplus/condor setup (due to port restrictions)",
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
        "--voms",
        default=None,
        type=str,
        help="Path to voms proxy, accessible to worker nodes. By default a copy will be made to $HOME.",
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
