"""Compare two strings"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "DV Klopfenstein, PhD"

from os.path import basename
from collections import namedtuple

def get_readpair_files(readpairs):
    """Get the base csv filename for each readpair"""
    nto = namedtuple('ReadPairs', 'fcsv fread1 fread2')
    nts = []
    for frd1, frd2 in readpairs:
        # TODO3: make robust if this is NOT a readpair
        fcsv = get_outputfname(frd1, frd2)
        nts.append(nto(fcsv=fcsv, fread1=frd1, fread2=frd2))
    return nts

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
    return res[0] if set(res[1:]) == {'1', '2'} else None

def get_outputfname(filename1, filename2):
    """Get the output file name given a read pair"""
    # TODO3: make robust if this is NOT a readpair
    pt0 = is_readpair(filename1, filename2)
    fcsv = basename(filename1)
    # TODO4: Ensure this is really an extension
    ptz = fcsv.find('.')
    if ptz != -1:
        fcsv = fcsv[:ptz]
    return f'{fcsv[:pt0]}p{fcsv[pt0+1:]}.csv'



# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
