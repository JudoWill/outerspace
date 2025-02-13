"""Extraction of sequence from reads"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "RE Berman"

from sys import exit as sys_exit
from os.path import isdir
from os.path import isfile
from os.path import exists
from argparse import ArgumentParser

def get_args(args=None):
    """cli for entering the reads"""
    parser = ArgumentParser(
                    prog='grna_extraction',
                    description='Get protospacers and UMIs',
                    epilog='Created by ThreeBlindMice - See how they run code')
    parser.add_argument('filename_or_dir', nargs='*',
            help='zipped fastq files or a directory containing fastq files')

    args = parser.parse_args(args)
    if args.filename_or_dir:
        _chk_exists(args.filename_or_dir)
    else:
        print('EXITING: NEED FASTQ READS FILES TO READ')
        sys_exit(1)
    print(f'ARGS: {args}')
    return args

def _chk_exists(filename_or_dir):
    not_exists = []
    for name in filename_or_dir:
        if not exists(name):
            not_exists.append(name)
        elif isdir(name):
            raise RuntimeError('TIME TO IMPLEMENT READING DIRECTORIES')
        else:
            assert isfile(name)
    if not_exists:
        for name in not_exists:
            print(f'DOES NOT EXIST: {name}')
        print(f'EXITING: {len(not_exists)} FASTQ FILE/DIR NOT EXIST')
        sys_exit(1)

# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
