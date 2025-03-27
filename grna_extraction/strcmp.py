"""Compare two strings"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "DV Klopfenstein, PhD"

# TODO1: from collections import namedtuple


def get_readpairs(files):
    """Return the files, grouped by read pairs, which differ by '1' vs '2'"""
    if not files:
        return None
    files = sorted(files)
    readpairs = []
    # TODO1: other = []
    # TODO1: return files that were not readpairs
    # TODO1: nto = namedtuple('ReadPairs', 'readpairs other')
    last = files[0]
    for curr in files[1:]:
        # TODO2: Optimize when a read pair is found
        if is_readpair(last, curr):
            readpairs.append((last, curr))
        last = curr
    # TODO1: return nto(readpairs=readpairs, other=other)
    return readpairs

def is_readpair(filename1, filename2):
    """Is a readpair if there is one difference between two files & it is '1' vs '2'"""
    if len(filename1) != len(filename2):
        return None
    res = [(i, a, b) for i, (a, b) in enumerate(zip(filename1, filename2)) if a != b]
    if len(res) != 1:
        return None
    res = res[0]
    return res if set(res[1:]) == {'1', '2'} else None


# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
