import pytest


def test_base_processor_warnings():
    """Test that warnings are raised when HggBaseProcessor is
    instantiated without arguments.
    The warnings are due to the processor not being able to find
    weights for PhotonIDMVA and DiPhotonMVA, since no metaconditions
    are provided.
    """
    from higgs_dna.workflows.base import HggBaseProcessor

    with pytest.warns(UserWarning):
        HggBaseProcessor()
