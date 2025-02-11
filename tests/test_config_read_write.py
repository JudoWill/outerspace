#!/usr/bin/env python3
"""testing config for read-write only"""

from os.path import exists
from grna_extraction.config import Cfg
print(f'TTTTTTTT test_config({__name__})')


def test_config():
    """testing config"""
    filename_config = "testing.cfg"
    cfg = Cfg()
    cfg.write_file(filename_config)
    assert exists(filename_config), f'DOES NOT EXIST: {filename_config}'


    print("TEST PASSED")

if __name__ == '__main__':
    test_config()
