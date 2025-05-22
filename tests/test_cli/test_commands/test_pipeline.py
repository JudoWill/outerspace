"""Tests for the pipeline command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
from argparse import Namespace
from outerspace.cli.commands.pipeline import PipelineCommand

@pytest.mark.skip(reason="Not implemented")
def test_pipeline_initialization():
    """Test that pipeline command initializes correctly"""
    args = Namespace(
        command='pipeline',
        config_filename=None,
        read1_filename=None,
        read2_filename=None,
        output_filename=None,
        fastqfiles=None,
        outdir=None
    )
    cmd = PipelineCommand(args)
    assert cmd.args == args

def test_pipeline_not_implemented():
    """Test that pipeline command is not yet implemented"""
    args = Namespace(
        command='pipeline',
        config_filename=None,
        read1_filename=None,
        read2_filename=None,
        output_filename=None,
        fastqfiles=None,
        outdir=None
    )
    cmd = PipelineCommand(args)
    with pytest.raises(NotImplementedError):
        cmd.run() 