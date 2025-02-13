#!/usr/bin/env python3
"""testing cli"""

from logging import debug
from logging import DEBUG
from logging import basicConfig
basicConfig (level = DEBUG)
from pytest import raises
from grna_extraction.cli import get_args

def test_cli():
    """testing no arguments"""
    with raises(SystemExit) as excinfo:
        args = get_args([])
        print(args)
    ## testing for error code 1 "systemexit 1" = missing read files
    ## to see if there are NO args -- we want an error here
    ## args with no FASTQ files should result in an error
    assert excinfo.value.code == 1
    print("TEST PASSED")

if __name__ == '__main__':
    test_cli()
