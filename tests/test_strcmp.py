#!/usr/bin/env python3
"""Test identifying files that are read pairs"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "DV Klopfenstein, PhD"

from itertools import permutations
from grna_extraction.strcmp import is_readpair
from grna_extraction.strcmp import get_readpairs


def test_strcmp():
    """Test finding read pairs"""
    files = [
        "a_R1_z.fq",
        "a_R2_z.fq",
    ]
    _prt_sorted(files)

    # Test a list of files that has one read pair
    assert is_readpair(*files)
    act = get_readpairs(files)
    exp = [tuple(files)]
    assert act == exp

    # Add files that are not read pairs
    files.append("a_RX_z.fq")
    files.append("apa")
    _prt_sorted(files)

    _run_readpairs(files)
    act = get_readpairs(files)
    assert act == exp
    print(f'ACT: {act}')
    print(f'EXP: {exp}')

    # Add files that are read pairs
    r12 = [
        "a_r1_z.fq",
        "a_r2_z.fq",
    ]
    files.extend(r12)
    exp.append(tuple(r12))
    _prt_sorted(files)


    act = get_readpairs(files)
    assert act == exp
    print(f'ACT: {act}')
    print(f'EXP: {exp}')
    #rps1 = [files[0], files[1]]
    #assert get_readpairs(files) == [tuple(files[1:])], 't3'


def _run_readpairs(files):
    for rd1, rd2 in permutations(files, r=2):
        is_rp = is_readpair(rd1, rd2)
        if {rd1, rd2} == {'a_R1_z.fq', 'a_R2_z.fq'}:
            assert is_rp is not None
            print(is_rp, rd1, rd2)
        else:
            assert is_rp is None, f'EXPECTED NOT READ PAIR: {rd1} {rd2}'

def _prt_sorted(files):
    print('')
    len_files = len(files)
    for idx, fin in enumerate(sorted(files)):
        print(f'{idx:2} {fin:9} files[{len_files}]')
    print('')


if __name__ == '__main__':
    test_strcmp()

# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
