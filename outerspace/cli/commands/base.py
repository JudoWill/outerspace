"""Base class for all CLI commands"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

from argparse import ArgumentParser
from sys import exit as sys_exit

class BaseCommand:
    """Base class for all commands"""
    def __init__(self, args=None):
        self.args = args

    def _init_parser(self, subparsers):
        """Initialize command-specific argument parser"""
        raise NotImplementedError("Each command must implement _init_parser")

    def run(self):
        """Execute the command"""
        raise NotImplementedError("Each command must implement run")


    def _chk_exists(self, filenames):
        """Check that files exist"""
        from os.path import exists, isdir, isfile
        not_exists = []
        for name in filenames:
            if not exists(name):
                not_exists.append(name)
            elif isdir(name):
                raise RuntimeError(f'Directory not supported: {name}')
            else:
                assert isfile(name)
        if not_exists:
            for name in not_exists:
                print(f'DOES NOT EXIST: {name}')
            print(f'EXITING: {len(not_exists)} FILE(S) NOT EXIST')
            sys_exit(1) 