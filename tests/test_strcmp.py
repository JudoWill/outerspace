#!/usr/bin/env python3
"""Test identifying files that are read pairs"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "DV Klopfenstein, PhD"

from itertools import permutations
from outerspace.strcmp import is_readpair
from outerspace.strcmp import get_readpairs
from outerspace.strcmp import get_outputfname
from outerspace.strcmp import get_readpair_files
from outerspace.strcmp import FastqSort

SEP = f"\n{'-'*80}\n"

def test_strcmp(pat=r'^(.*)(_R(1|2)_)(.*)$'):
    """Test finding read pairs"""
    print(f'{SEP}TEST 1) Two paired read files only')
    files = [
        "a_R1_z.fq.gz",
        "a_R2_z.fq.gz",
    ]
    _prt_sorted(files)

    # ------------------------------------------------------------
    # Test a list of files that has one read pair
    assert is_readpair(*files)
    act = get_readpairs(files, pat)
    exp = [tuple(files)]
    print(f'ACT: {act}')
    print(f'EXP: {exp}')
    assert act == exp
    fcsv = get_outputfname(*act[0])
    print(f'CSV({fcsv}) from readpair {files[0]} {files[1]}')
    assert fcsv == 'a_Rp_z.csv', fcsv

    # ------------------------------------------------------------
    # Add files that are not read pairs
    print(f'{SEP}TEST 2) Add files unrelated to reads')
    files.append("a_R1_32.fq")
    files.append("a_R1_34.fq")
    files.append("apa")
    _prt_sorted(files)

    #_run_readpairs(files)
    act = get_readpairs(files, pat)
    ntreadpairs = get_readpair_files(act)
    _prt_result(ntreadpairs)
    assert act == exp
    print(f'ACT: {act}')
    print(f'EXP: {exp}')

    # ------------------------------------------------------------
    # Add files that are read pairs
    print(f'{SEP}TEST 3) Add another read pair')

    r12 = [
        "a_r1_z.fq",
        "a_r2_z.fq",
    ]
    files.extend(r12)

    files.append("a_R2_32.fq")
    files.append("a_R2_34.fq")
    _prt_sorted(files)


    exp = [
        #"apa",             #  0
        ("a_r1_z.fq",       #  1
         "a_r2_z.fq"),      #  2
        ("a_R1_32.fq",      #  3
         "a_R2_32.fq"),     #  4
        ("a_R1_34.fq",      #  5
         "a_R2_34.fq"),     #  6
        ("a_R1_z.fq.gz",    #  7
         "a_R2_z.fq.gz"),   #  8
    ]
    act = get_readpairs(files, pat)
    assert act == exp, f'\nACT: {act}\nEXP: {exp}'
    print(f'ACT: {act}')
    print(f'EXP: {exp}')
    ntreadpairs = get_readpair_files(act)
    _prt_result(ntreadpairs)

def _run_readpairs(files):
    # TODO: Add expected results param
    mrks = set()
    for rd1, rd2 in permutations(files, r=2):
        is_rp = is_readpair(rd1, rd2)
        if {rd1, rd2} == {'a_R1_z.fq.gz', 'a_R2_z.fq.gz'}:
            assert is_rp is not None
            msg = f'FOUND READ PAIR: {is_rp}, {rd1}, {rd2}'
            if msg not in mrks:
                print(msg)
                mrks.add(msg)
        else:
            assert is_rp is None, f'EXPECTED NOT READ PAIR: {rd1} {rd2}'

def _prt_result(ntreadpairs):
    for ntd in ntreadpairs:
        print(ntd)

def _prt_sorted(files):
    print('')
    len_files = len(files)
    sorted_files = FastqSort().get_sorted(files)
    for idx, fin in enumerate(sorted_files):
        print(f'{idx:2} {fin:12} files[{len_files}]')
        #fintxt = f'"{fin}",'
        #print(f'        {fintxt:15}   # {idx:2}')
    print('')


if __name__ == '__main__':
    test_strcmp()

# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
