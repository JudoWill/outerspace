"""Tests for the main CLI functionality"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
from argparse import Namespace
from outerspace.cli.main import Cli


def test_cli_help(capsys):
    """Test that CLI shows help when no command is provided"""
    cli = Cli([])
    cli.run()
    captured = capsys.readouterr()
    assert "outerspace [-h] {findseq,collapse,count,gini,visualize,pipeline}" in captured.out

