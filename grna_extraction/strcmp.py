"""Compare two strings"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "DV Klopfenstein, PhD"

from os.path import basename
from operator import itemgetter
from collections import namedtuple
import re

def get_readpair_files(readpairs):
    """Get the base csv filename for each readpair"""
    nto = namedtuple('ReadPairs', 'fcsv fread1 fread2')
    nts = []
    for frd1, frd2 in readpairs:
        # TODO3: make robust if this is NOT a readpair
        fcsv = get_outputfname(frd1, frd2)
        nts.append(nto(fcsv=fcsv, fread1=frd1, fread2=frd2))
    return nts

def get_readpairs(files, pat=r'^(.*)(_R(1|2)_)(.*)$', seq=None):
    """Return the files, grouped by read pairs, which differ by '1' vs '2'"""
    obj = FastqSort(pat, seq)
    return obj.get_readpairs(files)

def is_readpair(filename1, filename2):
    """Is a readpair if there is one difference between two files & it is '1' vs '2'"""
    if len(filename1) != len(filename2):
        return None
    res = [(i, (a, b)) for i, (a, b) in enumerate(zip(filename1, filename2)) if a != b]
    if len(res) != 1:
        return None
    idx, set12 = res[0]
    return idx if set(set12) == {'1', '2'} else None

def get_outputfname(filename1, filename2):
    """Get the output file name given a read pair"""
    # TODO: make robust if this is NOT a readpair
    pt0 = is_readpair(filename1, filename2)
    fcsv = basename(filename1)
    # TODO: Ensure this is really an extension
    ptz = fcsv.find('.')
    if ptz != -1:
        fcsv = fcsv[:ptz]
    return f'{fcsv[:pt0]}p{fcsv[pt0+1:]}.csv'


class FastqSort:
    """Sort fastq files such that read1 and read2 files are next to each other"""

    # TODO: Make this work for the world, not just Dampier
    default_dampier_pat = r'^(.*)(_R(1|2)_)(.*)$'
    default_dampier_seq = [0, 3, 1]
    default_sort_nomatch = ()


    def __init__(self, pat=r'^(.*)(_R(1|2)_)(.*)$', seq=None):
        # TODO: Make more robust for various patterns and sequences
        self.pat = pat
        self.seq = self._init_seq(pat, seq)
        self.cmp = re.compile(pat)
        self.getter = itemgetter(*self.seq)

    def get_sorted(self, files):
        """Get files sorted to match paired reads together"""
        return sorted(files, key=self._key)

    def get_readpairs(self, files):
        """Return the files, grouped by read pairs, which differ by '1' vs '2'"""
        if not files:
            return None
        files = self.get_sorted(files)
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

    def _key(self, item):
        if mtch := self.cmp.search(item):
            return self.getter(mtch.groups())
        return self.default_sort_nomatch

    def _init_seq(self, pat, seq):
        if seq is None:
            if pat == self.default_dampier_pat:
                return self.default_dampier_seq
            raise RuntimeError("Please provide correct sequence for regex pat({pat})")
        # TODO: Check that seq is a list (or tuple)
        return seq


# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
