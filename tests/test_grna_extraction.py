#!/usr/bin/env python3
"""testing top level read pattern matching"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "??"


from os import remove
from os.path import exists
from os.path import join
from grna_extraction.grna_extraction import TopLevel
from grna_extraction.read_fastq import ReadPairedFastq
from tests.pkgtest.utils import mk_outdir
from tests.pkgtest.utils import get_filename
print(f'TTTTTTTT test_search({__name__})')

#Testing the file that runs all the modules
def test_grna_extraction():
    """testing config"""
    # TODO: replace with test data
    # Output path and out file name
    # Setting a variable to what will be the name
    fcsv = join(mk_outdir('outdir'), 'out.csv')
    # TODO: fq1 needs to be copied over again- messed it up
    # Input path
    # fq1 = "/data/share/nonn-lab/rachel-test-crispr/reads/409-4_S1_L001_R1_001.fastq.gz"
    # fq2 = "/data/share/nonn-lab/rachel-test-crispr/reads/409-4_S1_L001_R2_001.fastq.gz"
    fq1 = get_filename("tests/data/409-4_S1_L002_R1_001.fastq.gz")
    fq2 = get_filename("tests/data/409-4_S1_L002_R2_001.fastq.gz")
    assert exists(fq1)
    assert exists(fq2)
    cfg_filename = get_filename("tests/configs/grnaquery.cfg")
    assert exists(cfg_filename)
    
    obj = TopLevel(cfg_filename) 
    # if file exists remove it so that new file can take its place
    if exists(fcsv):
        remove(fcsv)
    assert not exists(fcsv)
    obj.run(fcsv, fq1, fq2, do_break=True)
    # Creating configuration object
    # Reads configuration files and loads that info into dictionary and delivers that back to you

    assert obj.search.srch.regxlist == []
    assert obj.search.srch1.regxlist == [
        '(?P<UMI_5prime>.{8})(?:CTTGGCTTTATATATCTTGTGG){s<=4}',
        '(?:TATCTTGTGGAAAGGACGAAACACC){s<=4}(?P<protospacer>.{19,21})(?P<downstreamof_protospacer>GTTTAAGTACTCTGTGCTGGAAACAG){s<=4}'
        ]
    assert obj.search.srch2.regxlist ==  ['(?P<UMI_3prime>.{8})(?:TTCCACACCCTAACTGACACAC){s<=4}']

    # Is there anything in list
    assert obj.search.srch1.cmps
    # THe compiler found capture names in object
    assert obj.search.srch1.capturednames
    print(f'NAMES1: {obj.search.srch1.names}')

    assert exists(fq1)
    assert exists(fq2)

    # We are writing out to the csv file, and we are using both reads
    # Setting output and giving files that contain fastq giving it to reader so that it iterates through them one by one
    assert exists(fcsv)
    #TODO: Add a check for the printed output file to make sure its correctly printing stuff out
    print("TEST PASSED")



if __name__ == '__main__':
    test_grna_extraction()


# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
