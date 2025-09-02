#!/usr/bin/env python3
"""Test uppercase and reverse complement a nuceotide pattern"""

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

def reverse_complement(seq):
    """Reverse complement a DNA sequence"""
    complement = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "a": "t",
        "c": "g",
        "g": "c",
        "t": "a",
    }
    return "".join(complement.get(base, base) for base in reversed(seq))


def test_revpat():
    """Test uppercase and reverse complement a nuceotide pattern"""
    back_umi_forward = "gtgtgtcagttagggtgtggaa".upper()
    assert back_umi_forward == "GTGTGTCAGTTAGGGTGTGGAA"
    back_umi_rc = reverse_complement(back_umi_forward)
    assert back_umi_rc == "TTCCACACCCTAACTGACACAC"
    pat = f"(?P<UMI>.{{8}})(?:{back_umi_rc}){{s<=4}}"
    assert pat == "(?P<UMI>.{8})(?:TTCCACACCCTAACTGACACAC){s<=4}"
    print(pat)


if __name__ == "__main__":
    test_revpat()

    
# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.