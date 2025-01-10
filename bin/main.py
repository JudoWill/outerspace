#!/usr/bin/env python3
"""Extraction of sequence from reads"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "SC Barrera"

from os.path import join
#from grna_extraction.extraction_attempt import process_paired_read_file

from Bio.Seq import reverse_complement

import regex


def main():
    print('running_main')
    csv = '409-4.csv'

    dir_crispr = '../../../nonn-lab/rachel-test-crispr'

    dir_read = join(dir_crispr, 'reads/')
    path1 = join(dir_read, '409-4_S1_L001_R1_001.fastq.gz')
    path2 = join(dir_read, '409-4_S1_L001_R2_001.fastq.gz')

    forward_umi_reg = regex.compile('(?P<UMI>.{8})(?:CTTGGCTTTATATATCTTGTGG){s<=4}', flags=regex.BESTMATCH) 
    protospacer_reg = regex.compile('(?:TATCTTGTGGAAAGGACGAAACACC){s<=4}(?P<protospacer>.{19,21})(?:GTTTAAGTACTCTGTGCTGGAAACAG){s<=4}', flags=regex.BESTMATCH)
    back_umi_forward = 'gtgtgtcagttagggtgtggaa'.upper()
    back_umi_rc = reverse_complement(back_umi_forward)
    reverse_umi_reg = regex.compile(f'(?P<UMI>.{{8}})(?:{back_umi_rc}){{s<=4}}', flags=regex.BESTMATCH) 

    print(f'path1 = {path1}')
    print(f'path2 = {path2}')
    #process_paired_read_file(csv, path1, path2, forward_umi_reg, protospacer_reg, reverse_umi_reg)

if __name__ == '__main__':
    main()

# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
