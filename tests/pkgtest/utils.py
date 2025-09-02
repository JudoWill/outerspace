"""Utilities for test"""

# from os.path import exists
from os import makedirs
from os.path import join
from os.path import dirname
from os.path import normpath

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "DVK"

def get_filename(fname):
    """making absolute path regardless of what machine you are on"""
    return normpath(join(dirname(__file__), "../..", fname))


def mk_outdir(outdir):
    """Making output directory using absolute path name"""
    # getsfilename before path exists
    abs_outdir = get_filename(outdir)
    # makes directory, no error is raised if it already exists- we want this
    makedirs(abs_outdir, exist_ok=True)
    print(f"makedirs({abs_outdir})")
    return abs_outdir


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.