from higgs_dna.samples.fetch import get_dataset_dict
import pytest
from subprocess import getstatusoutput


@pytest.mark.skipif(
    getstatusoutput("dasgoclient --help")[0]
    or getstatusoutput("voms-proxy-info -e -p")[0],
    reason="dasgoclient or voms not installed",
)
def test_get_dataset_dict():
    lst = [
        (
            "DoubleEG-Run2017B",
            "/DoubleEG/Run2017B-UL2017_MiniAODv1_NanoAODv2-v1/NANOAOD",
        ),
        (
            "DYJets-M50",
            "/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv2-106X_mc2017_realistic_v8-v1/NANOAODSIM",
        ),
    ]

    samples = get_dataset_dict(lst, "root://cms-xrd-global.cern.ch/", "prod/global")

    assert "DoubleEG-Run2017B" in samples
    assert "DYJets-M50" in samples
