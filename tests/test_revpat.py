#!/usr/bin/env python3
"""Test uppercase and reverse complement a nuceotide pattern"""

import Bio
from Bio.Seq import reverse_complement


def test_revpat():
    """Test uppercase and reverse complement a nuceotide pattern"""
    back_umi_forward = 'gtgtgtcagttagggtgtggaa'.upper()
    assert back_umi_forward == 'GTGTGTCAGTTAGGGTGTGGAA'
    back_umi_rc = reverse_complement(back_umi_forward)
    assert back_umi_rc == 'TTCCACACCCTAACTGACACAC'
    pat = f'(?P<UMI>.{{8}})(?:{back_umi_rc}){{s<=4}}'
    assert pat == '(?P<UMI>.{8})(?:TTCCACACCCTAACTGACACAC){s<=4}'
    print(pat)


if __name__ == '__main__':
    test_revpat()
