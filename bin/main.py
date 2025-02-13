#!/usr/bin/env python3
"""Extraction of sequence from reads"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "SC Barrera"

from sys import exit as sys_exit
from os.path import join
from os.path import isdir
from os.path import isfile
from os.path import exists
from argparse import ArgumentParser

from Bio.Seq import reverse_complement

import regex

from original.extraction_attempt import process_paired_read_file

def main():
    print('running_main')
    csv = '409-4.csv'

    ####dir_crispr = '../../../nonn-lab/rachel-test-crispr'

    ####dir_read = join(dir_crispr, 'reads/')
    args = _get_args()
    ####path1 = join(dir_read, '409-4_S1_L001_R1_001.fastq.gz')
    ####path2 = join(dir_read, '409-4_S1_L001_R2_001.fastq.gz')
    path1 = args.filename_or_dir[0]
    path2 = args.filename_or_dir[1]

    forward_umi_reg = regex.compile('(?P<UMI>.{8})(?:CTTGGCTTTATATATCTTGTGG){s<=4}', flags=regex.BESTMATCH) 
    protospacer_reg = regex.compile('(?:TATCTTGTGGAAAGGACGAAACACC){s<=4}(?P<protospacer>.{19,21})(?:GTTTAAGTACTCTGTGCTGGAAACAG){s<=4}', flags=regex.BESTMATCH)
    back_umi_forward = 'gtgtgtcagttagggtgtggaa'.upper()
    back_umi_rc = reverse_complement(back_umi_forward)
    reverse_umi_reg = regex.compile(f'(?P<UMI>.{{8}})(?:{back_umi_rc}){{s<=4}}', flags=regex.BESTMATCH) 

    print(f'path1 = {path1}')
    print(f'path2 = {path2}')
    print(f'forward_umi_reg = {forward_umi_reg.pattern}')
    print(f'protospacer_reg = {protospacer_reg.pattern}')
    print(f'reverse_umi_reg = {reverse_umi_reg.pattern}')
    
    process_paired_read_file(csv, path1, path2, forward_umi_reg, protospacer_reg, reverse_umi_reg)

def _get_args():
    parser = ArgumentParser(
                    prog='grna_extraction',
                    description='Get protospacers and UMIs',
                    epilog='Created by ThreeBlindMice - See how they run code')
    parser.add_argument('filename_or_dir', nargs='*',
            help='zipped fastq files or a directory containing fastq files')
    
    args = parser.parse_args()
    if args.filename_or_dir:
        _chk_exists(args.filename_or_dir)
    else:
        print(f'EXITING: NEED FASTQ READS FILES TO READ')
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
        
if __name__ == '__main__':
    main()

# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
