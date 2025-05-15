"""Parse reads, find user-given patterns, save sequence matches to csv."""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "SC Barrera"


from os.path import exists
from outerspace.config import Cfg
from outerspace.top_search import TopSearch
from outerspace.read_fastq import ReadPairedFastq


# This function runs through all package modules to parse the reads, compiles patterns using regx, finds requested patterns within sequences, then writes to a csv file
# Give cfg file, fastq file1 fastq file2, output filename; if not given output to stdout

def run(file_pathname_cfg, file_pathname_fastq1, file_pathname_fastq2=None, filename_csv=None):
    """Parse reads, find user-given patterns, save sequence matches to csv."""
    obj = TopLevel(file_pathname_cfg)
    obj.run(file_pathname_fastq1, file_pathname_fastq2, filename_csv)

class TopLevel:
    """Parse reads, find user-given patterns, save sequence matches to csv."""
    # TODO: pick a built in python exception and informative message
    # Making an object that contains all pieces
    def __init__(self, file_pathname_cfg):

        # Output path and out file name, Setting a variable to what will be the name
        assert exists(file_pathname_cfg), file_pathname_cfg

        # Creating configuration object
        # Reads configuration files and loads that info into dictionary and delivers that back to you
        self.cfg = Cfg(file_pathname_cfg)
        self.doc = self.cfg.read_file()

        #
        self.search = TopSearch(self.doc)

        # TODO: Make this a nicer check with pleasant messages to researcher
        # Is there anything in list
        assert self.search.srch1.cmps
        # THe compiler found capture names in object
        assert self.search.srch1.capturednames
        # creating an object- giving it search (not yet reading), ReadPairedFastq is from read fastq file
        self.reader = ReadPairedFastq(self.search)

    def run(self, filename_csv, file_pathname_fastq1, file_pathname_fastq2, do_break=False):
        """Parse reads, find user-given patterns, save sequence matches to csv."""
        # TODO: Make this nicer - used as a communication tool
        assert exists(file_pathname_fastq1), file_pathname_fastq1

        # We are writing out to the csv file, and we are using both reads
        # Setting output and giving files that contain fastq giving it to reader so that it iterates through them one by one
        self.reader.process_paired_read_file(filename_csv, file_pathname_fastq1,file_pathname_fastq2, do_break)
        #TODO: Add a check for the printed output file to make sure its correctly printing stuff out


# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
