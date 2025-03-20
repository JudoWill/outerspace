"""Extraction of sequence from reads"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "RE Berman"

from sys import argv
from sys import exit as sys_exit
from os.path import isdir
from os.path import isfile
from os.path import exists
from argparse import ArgumentParser
from grna_extraction.grna_extraction import TopLevel

class Cli:
    """The command line interface argument requirements"""
    def __init__(self, args=None):
        self.parser = self._init_parser()
        print(f'ARGV: {argv}')
        self.args = self.parser.parse_args(args)
        print(f'ARGS: {self.args}')
        self.top = self._init_top()
    
    def run(self):
        """running if there are reads"""
        # there are four conditions 2 reads, no reads, read 1 and not read 2, read 2 and not read1
        if self.args.read1_filename is not None and self.args.read2_filename is not None:
            self._runpairedreads()
        elif self.args.read1_filename is not None and self.args.read2_filename is None:
            self._runsingleread()
        else:
            self._runnone()
 
    def _runpairedreads(self):
        """running if there are paired reads files"""
        # TODO: check that reads files exist and that the output file was given, if it wasnt given we should store it
        self.top.run(self.args.output_filename, self.args.read1_filename, self.args.read2_filename)
        # raise NotImplementedError("TIME TO IMPLEMENT RUNNING THE PAIRED READS")
    
    def _runsingleread(self):
        raise NotImplementedError("TIME TO IMPLEMENT RUNNING SINGLE READ")
    
    def _runnone(self):
        raise NotImplementedError("TIME TO EXIT - NO READ FILES")
        
    @staticmethod
    def _init_parser():
        """cli for entering the reads"""
        parser = ArgumentParser(
                        prog='grna_extraction',
                        description='Get protospacers and UMIs',
                        epilog='Created by ThreeBlindMice - See how they run code')
        parser.add_argument('config_filename', nargs='?',
            help='Configuration file holds user input search patterns ')
        parser.add_argument('-1', '--read1_filename',
            help='zipped fastq file for read 1, or a single read')
        parser.add_argument('-2', '--read2_filename',
            help='zipped fastq file for read 2, or a single read')
        parser.add_argument('-o', '--output_filename',
            help='captured read file name output csv')
        return parser

    @staticmethod
    def _chk_exists(config_filename):
        """checking that reads exist"""
        not_exists = []
        for name in config_filename:
            if not exists(name):
                not_exists.append(name)
            elif isdir(name):
                raise RuntimeError('TIME TO IMPLEMENT READING DIRECTORIES')
            else:
                assert isfile(name)
        if not_exists:
            for name in not_exists:
                print(f'DOES NOT EXIST: {name}')
            print(f'EXITING: {len(not_exists)} FASTQ FILE/DIR NOT EXIST')
            sys_exit(1)

    def _init_top(self):
        """checking that config exists"""
        filenamecfg = self.args.config_filename
        if filenamecfg is None:
            # TODO: give researcher explicit instructions how to generate a config file - give hints or examples
            print("Provide configuration file name")
            return None
        if exists(filenamecfg):
            return TopLevel(filenamecfg)
        print(f'Configuration file does not exist: {filenamecfg}')
        return None
    

# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
