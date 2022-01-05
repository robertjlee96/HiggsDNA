import pytest
from higgs_dna.workflows import workflows

def test_processor_behavior():
    fake_dict = {"meta": None}

    # Create instance of base processor with dictionary which does not correspond to any metacondition
    # Will raise a KeyError
    with pytest.raises(KeyError):
        workflows["dystudies"](fake_dict)