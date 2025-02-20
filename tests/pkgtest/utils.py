"""Utilities for test"""

# from os.path import exists
from os.path import join
from os.path import dirname
from os.path import normpath

def get_filename(fname):
    """ making absolute path regardless of what machine you are on"""
    return normpath(join(dirname(__file__), '../..',fname))
