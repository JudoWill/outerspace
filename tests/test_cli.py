#!/usr/bin/env python3
"""testing cli"""

from grna_extraction.cli import get_args

def test_cli():
    """testing cli"""
    args = get_args()
    print(args)

    print("TEST PASSED")

if __name__ == '__main__':
    test_cli()
