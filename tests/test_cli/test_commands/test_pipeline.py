"""Tests for the pipeline command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
from outerspace.cli.main import Cli

@pytest.mark.skip(reason="Not implemented")
def test_pipeline_initialization():
    """Test that pipeline command initializes correctly"""
    args = [
        'pipeline'
    ]
    cli = Cli(args)
    assert cli.args.command == 'pipeline'

def test_pipeline_not_implemented():
    """Test that pipeline command is not yet implemented"""
    args = [
        'pipeline'
    ]
    cli = Cli(args)
    with pytest.raises(NotImplementedError):
        cli.run() 