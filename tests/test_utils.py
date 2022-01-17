import pytest


def test_logger():
    """
    Test main aspects of setup_logger function
    """
    from higgs_dna.utils.logger_utils import setup_logger

    null_info_level = "NOT_ALLOWED"
    with pytest.raises(ValueError):
        setup_logger(null_info_level)


def test_parser():
    """
    Test the parser

    1. Check that passing no arguments raises SystemExit
    2. Check that passing only part of the required arguments raises SystemExit
    3. Check that passing a non-existing argument raises SystemExit
    4. Check that fields are correctly parsed
    """
    from higgs_dna.utils.runner_utils import get_main_parser

    parser = get_main_parser()

    with pytest.raises(SystemExit):
        parser.parse_args([])

    with pytest.raises(SystemExit):
        parser.parse_args(["--wf", "dystudies"])

    with pytest.raises(SystemExit):
        parser.parse_args(["--fake-arg"])

    args = parser.parse_args(["--wf", "dystudies", "--meta", "Era2017_legacy_xgb_v1"])
    assert args.workflow == "dystudies"
    assert args.metaconditions == "Era2017_legacy_xgb_v1"
