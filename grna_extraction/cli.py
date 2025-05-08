"""Extraction of sequence from reads"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "RE Berman"

# from sys import argv
from sys import exit as sys_exit
from os.path import isdir
from os.path import isfile
from os.path import exists
from os.path import join
from argparse import ArgumentParser
from timeit import default_timer
from datetime import timedelta
from grna_extraction.grna_extraction import TopLevel
from grna_extraction.strcmp import get_readpair_files
from grna_extraction.strcmp import get_readpairs

def main(args=None):
    """Main runs outerspace- run from setup.py console scripts """
    cli=Cli(args)
    cli.run()

class Cli:
    """The command line interface argument requirements"""
    def __init__(self, args=None):
        self.parser = self._init_parser()
        # print(f'ARGV: {argv}')
        self.args = self.parser.parse_args(args)
        self._prt_args()
        self.top = self._init_top()
    
    def run(self):
        """running if there are reads"""
        if self.args.config_filename is None:
            print('please provide a config filename')
            sys_exit(0)
        if self.args.read1_filename is not None and \
           self.args.read2_filename is not None and \
           self.args.output_filename is not None:
            self._run_reads_one(self.args.read1_filename,self.args.read2_filename)
        if self.args.fastqfiles and self.args.outdir:
            self._run_pairedreads_all(self.args.fastqfiles,self.args.outdir)
    
    def _toprun(self, outputfile, read1file, read2file):
        """putting in a clause that will make this work even if theres an error"""
        try:
            self.top.run(outputfile, read1file, read2file)
        except ValueError as err:
            print(f'Failed analyzing read pair: {read1file} {read2file}')
            
    def _run_pairedreads_all(self, fastqfiles, outdir):
        """running if there is a directory with all the paired reads you intend to run"""
        # testing to see if outdir exists
        if not exists(outdir) or not isdir(outdir):
            print('Output directory either does not exist or is not a directory:({outdir})')
            sys_exit(0)
        if not fastqfiles:
            print('No fastq files were found with --fastqfiles {fastqfiles}')
            sys_exit(0)
        tic = default_timer()
        readpairs = get_readpairs(fastqfiles)
        # TODO: check that we have nts (to check that we have readpairs)
        nts = get_readpair_files(readpairs)
        len_nts = len(nts)
        for idx,ntd in enumerate(nts,1):
            hms = timedelta(seconds = default_timer() - tic)
            print(f'hours:mins:secs {hms} -- {idx:4} of {len_nts} readpair {ntd}')
            output = join(outdir, ntd.fcsv)
            self._toprun(output, ntd.fread1, ntd.fread2)

    def _run_reads_one(self, read1_filename,read2_filename):    
        # there are four conditions 2 reads, no reads, read 1 and not read 2, read 2 and not read1
        if read1_filename is not None and read2_filename is not None:
            self._runpairedreads()
        elif read1_filename is not None and read2_filename is None:
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
        print("TODO: make exit friendly")
        
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
        parser.add_argument('--fastqfiles', nargs='*',
            help='directory containing paired fastq read files')
            # TODO: add option to specify single reads vs paired reads
        parser.add_argument('--outdir',
            help='output directory for processed files')
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
            #TODO:consider if want sys_exit(0); report to log file
            sys_exit(1)

    def _init_top(self):
        """checking that config exists"""
        filenamecfg = self.args.config_filename
        if filenamecfg is None:
            # TODO: give researcher explicit instructions how to generate a config file - give hints or examples
#            print("Provide configuration file name")
            return None
        if exists(filenamecfg):
            return TopLevel(filenamecfg)
        print(f'Configuration file does not exist: {filenamecfg}')
        return None
    
    def _prt_args(self):
        """shortening the output of args to summarize input r1 r2 files, output file or dir, and/or the amount of files found in the input dir"""
        txt = []
        for key,val in vars(self.args).items():
            if key != 'fastqfiles':
                txt.append(f'{key}={val}')
            else:
                txt.append(f'{key}[{0 if val is None else len(val)}]')
        txt = ', '.join(txt)
        print(f'ARGS: {txt}')
    

# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
