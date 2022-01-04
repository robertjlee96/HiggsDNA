import pytest
from higgs_dna.utils.logger_utils import setup_logger


def test_logger():
    """
    Test main aspects of setup_logger function
    """
    null_info_level = "NOT_ALLOWED"
    with pytest.raises(ValueError):
        setup_logger(null_info_level)