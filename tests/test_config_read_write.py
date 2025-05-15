#!/usr/bin/env python3
"""testing config for read-write only"""

from os.path import exists
from os.path import join  
from os.path import dirname
from os.path import normpath
from outerspace.config import Cfg
print(f'TTTTTTTT test_config({__name__})')

def get_filename(fname):
    # making absolute path regardless of what machine on
    return normpath(join(dirname(__file__), '..',fname))

def test_config():
    """testing config"""
    filename_config = get_filename("example_file.cfg")
    cfg = Cfg(filename_config)
    cfg.write_file()
    doc = cfg.read_file()
    print(doc['regxlist'])
    assert exists(filename_config), f'DOES NOT EXIST: {filename_config}'


    print("TEST PASSED")

if __name__ == '__main__':
    test_config()
