#!/usr/bin/env python3
"""testing config"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "??"


from os import remove
from os.path import exists
from os.path import join
from grna_extraction.config import Cfg
from grna_extraction.search import Search
from grna_extraction.read_fastq import ReadPairedFastq
from tests.pkgtest.utils import mk_outdir
print(f'TTTTTTTT test_search({__name__})')


def test_read_fastq():
    """testing config"""
    # TODO: replace with test data
    # Output path and out file name
    fcsv = join(mk_outdir('outdir'), 'out.csv') 
    # Input path
    fq1 = "/data/share/nonn-lab/rachel-test-crispr/reads/409-4_S1_L001_R1_001.fastq.gz"
    fq2 = "/data/share/nonn-lab/rachel-test-crispr/reads/409-4_S1_L001_R2_001.fastq.gz"
    
    cfg = Cfg()
    doc = cfg.get_doc_default()
    regxlist1 = doc['regxlist1']
    search = Search(regxlist1)
    
    print(f'NAMES1: {search.capture_names()}')
    assert search.capture_names()
    assert search.cmps
    reader = ReadPairedFastq(search)

    assert exists(fq1)
    assert exists(fq2)
    # if file exists remove it so that new file can take its place
    if exists(fcsv):
        remove(fcsv)
    assert not exists(fcsv)

    # We are writing out to the csv file, and we are using both reads
    reader.process_paired_read_file(fcsv, fq1,fq2)
    assert exists(fcsv)
    #TODO: Add a check for the printed output file to make sure its correctly printing stuff out
    print("TEST PASSED")

    

if __name__ == '__main__':
    test_read_fastq()





 # Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
