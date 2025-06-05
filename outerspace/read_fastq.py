"""Extraction of sequence from reads"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "SC Barrera"

from sys import exit as sys_exit
import gzip
from tqdm import tqdm
import pyfastx
import csv
from collections import namedtuple

# Define Read as a namedtuple at module level
Read = namedtuple('Read', ['id', 'seq', 'qual'])

class ReadPairedFastq:
    """Reading in & parsing the fastq paired reads"""
    def __init__(self, topsearch):
        self.topsearch = topsearch

    def _extract_from_paired_reads(self, file1, file2, do_break=False):
        """Actually reading & parsing the fastq paired reads from given files"""
        total = 0
        missed = 0

        #enumerate starting at 1
        topsearch = self.topsearch
        for total, (read1, read2) in enumerate(tqdm(self._iterate_readpairs(file1, file2)),1):
            # topsearch is getting the capture from read1 and read2, calling it matches
            matches = topsearch.get_capture_from_readpair(read1,read2)

            if matches:
                yield matches
            else:
                missed +=1
            # if doing break and more than 10 reads then breaking
            if do_break and total > 10:
                break

    @staticmethod
    def _iterate_reads(path):
        "Iterate reads from a gzipped or regular fastq files"
        # pyfastx automatically handles gzipped files
        for read in pyfastx.Fastx(path, build_index=False):
            # Convert pyfastx read to a format compatible with existing code
            # pyfastx returns (name, sequence, quality) tuple
            yield Read(read[0], read[1], read[2])

    def _iterate_readpairs(self, path1, path2):
        for r1, r2 in zip(self._iterate_reads(path1), self._iterate_reads(path2)):
            yield r1, r2

    def process_paired_read_file(self, outpath, path1, path2, do_break=False):
        """Function takes headers and writes them"""
        fieldnames = self.topsearch.get_names_readpair()
        if not fieldnames:
            # TODO: want to make log file with various mesages for ppl
            # Exit if no capture names
            sys_exit(0)

        with open(outpath, mode='w') as handle:
            # Access the search object and get names
            writer = csv.DictWriter(handle, fieldnames)
            writer.writeheader()
            # dobreak is passing whatever the value is that is passed into function
            stream = self._extract_from_paired_reads(path1, path2, do_break)
            writer.writerows(stream)
            #TODO: Add logging statement of how many total reads processed and how may did not have pattern
