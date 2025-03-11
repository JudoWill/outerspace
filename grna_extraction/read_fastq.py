"""Extraction of sequence from reads"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "SC Barrera"

## ####from os.path import join
from sys import exit as sys_exit
import gzip
from tqdm import tqdm
from Bio import SeqIO

## ####from itertools import islice
## ####from Bio.Seq import reverse_complement
## ####
## ####import regex
import csv

class ReadPairedFastq:
    """Reading in & parsing the fastq paired reads"""
    def __init__(self, topsearch):
        self.topsearch = topsearch

    def _extract_from_paired_reads(self, file1, file2, do_break=False):
        """Actually reading & parsing the fastq paired reads from given files"""
        #def extract_from_paired_reads(read1_path, read2_path, forward_reg, proto_reg, reverse_reg):
        total = 0
        missed = 0

        #enumerate starting at 1
        topsearch = self.topsearch
        for total, (read1, read2) in enumerate(tqdm(self._iterate_readpairs(file1, file2)),1):
            # topsearch is getting the capture from read1 and read2, calling it matches
            matches = topsearch.get_capture_from_readpair(read1,read2)

            #function call, definition below for _process_matches
            # self._process_matches(1, read1, total, names1, matches1)
            if matches:
                yield matches
            else:
                missed +=1
            # if doing break and more than 10 reads then breaking
            if do_break and total > 10:
                break
        ####        foward_umi = get_capture_from_read(forward_reg, read1)
        ####        protospacer = get_capture_from_read(proto_reg, read1)
        ####        reverse_umi = get_capture_from_read(reverse_reg, read2)
        ####
        ####        if foward_umi and protospacer and reverse_umi:
        ####            yield {'forward_umi': foward_umi,
        ####                   'protospacer': protospacer,
        ####                   'reverse_umi': reverse_umi,
        ####                   'read_id': read1.id
        ####                  }
        ####        else:
        ####            missed += 1
        print(f'Total sequences: {total:,}\nMissed sequences: {missed}')

    @staticmethod
    def _process_matches(r1_2, read, readnum, names, matches):
        pass


    @staticmethod
    def _iterate_reads(path):
        "Iterate reads from a gzipped or regular fastq files"

        if path.endswith('.gz'):
            with gzip.open(path, mode='rt') as handle:
                for read in SeqIO.parse(handle, 'fastq'):
                    yield read
        else:
            with open(path) as handle:
                for read in SeqIO.parse(handle, 'fastq'):
                    yield read

    def _iterate_readpairs(self, path1, path2):
        for r1, r2 in zip(self._iterate_reads(path1), self._iterate_reads(path2)):
            yield r1, r2

## ####def get_capture_from_read(reg_exp, read):
## ####
## ####    sequence = str(read.seq)
## ####    result = reg_exp.findall(sequence)
## ####
## ####    if len(result) == 1:
## ####        return result[0]
## ####
    # Function takes headers  and writes them
    def process_paired_read_file(self, outpath, path1, path2):
        """Function takes headers  and writes them"""
        fieldnames = self.topsearch.get_objsearch('read1').get_csvheaders()
        print(f'CAPTURED FIELD NAMES: {fieldnames}')
        if not fieldnames:
            print(f'NO CAPTURE NAMES IN SEARCH PATTERNS')
            # TODO: want to make log file with various mesages for ppl
            # Exit if no capture names
            sys_exit(0)

        with open(outpath, mode='w') as handle:
            # Access the search object and get names
            writer = csv.DictWriter(handle, fieldnames)
            writer.writeheader()
            stream = self._extract_from_paired_reads(path1, path2)
            writer.writerows(stream)
            print(f' FILE WRITTEN: {outpath}')
            #TODO: Add print statement of how many total reads processed and how may did not have pattern
##
## ####forward_umi_reg = regex.compile('(?P<UMI>.{8})(?:CTTGGCTTTATATATCTTGTGG){s<=4}', flags=regex.BESTMATCH)
## ####protospacer_reg = regex.compile('(?:TATCTTGTGGAAAGGACGAAACACC){s<=4}(?P<protospacer>.{19,21})(?:GTTTAAGTACTCTGTGCTGGAAACAG){s<=4}',
## ####                                flags=regex.BESTMATCH)
## ####
## ####back_umi_forward = 'gtgtgtcagttagggtgtggaa'.upper()
## ####back_umi_rc = reverse_complement(back_umi_forward)
## ####
## ####reverse_umi_reg = regex.compile(f'(?P<UMI>.{{8}})(?:{back_umi_rc}){{s<=4}}', flags=regex.BESTMATCH)
## ####
## ##### Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
