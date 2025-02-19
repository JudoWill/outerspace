#!/usr/bin/env python3
"""testing config"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "??"


from os.path import exists
from grna_extraction.config import Cfg
from grna_extraction.search import Search
from grna_extraction.read_fastq import ReadPairedFastq
print(f'TTTTTTTT test_search({__name__})')


def test_read_fastq():
    """testing config"""
    # TODO: replace with test data
    fq1 = "/data/share/nonn-lab/rachel-test-crispr/reads/409-4_S1_L001_R1_001.fastq.gz"
    fq2 = "/data/share/nonn-lab/rachel-test-crispr/reads/409-4_S1_L001_R2_001.fastq.gz"
    cfg = Cfg()
    doc = cfg.get_doc_default()
    regxlist = doc['regxlist']
    search = Search(regxlist)
    reader = ReadPairedFastq(search)

    assert exists(fq1)
    assert exists(fq2)
    reader.run(fq1,fq2)
    print("TEST PASSED")
    

if __name__ == '__main__':
    test_read_fastq()





 # Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
