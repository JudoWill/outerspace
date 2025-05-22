"""Sequence extraction command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

from argparse import ArgumentParser
from os.path import join, exists, isdir
from timeit import default_timer
from datetime import timedelta
from sys import exit as sys_exit

from outerspace.cli.commands.base import BaseCommand
from outerspace.outerspace import TopLevel
from outerspace.strcmp import get_readpair_files, get_readpairs

class FindSeqCommand(BaseCommand):
    """Command for extracting sequences from FASTQ files"""
    def _init_parser(self, subparsers):
        """Initialize command-specific argument parser"""
        parser = subparsers.add_parser('findseq',
            help='Extract sequences from FASTQ files based on configuration patterns')
        parser.add_argument('config_filename', nargs='?',
            help='Configuration file with search patterns')
        parser.add_argument('-1', '--read1_filename',
            help='Zipped FASTQ file for read 1, or a single read')
        parser.add_argument('-2', '--read2_filename',
            help='Zipped FASTQ file for read 2, or a single read')
        parser.add_argument('-o', '--output_filename',
            help='Captured read file name output CSV')
        parser.add_argument('--fastqfiles', nargs='*',
            help='Directory containing paired FASTQ read files')
        parser.add_argument('--outdir',
            help='Output directory for processed files')
        return parser

    def run(self):
        """Run the findseq command"""
        if self.args.config_filename is None:
            print('Please provide a config filename')
            sys_exit(0)

        if self.args.read1_filename is not None and \
           self.args.read2_filename is not None and \
           self.args.output_filename is not None:
            self._run_reads_one(self.args.read1_filename, self.args.read2_filename)
        elif self.args.fastqfiles and self.args.outdir:
            self._run_pairedreads_all(self.args.fastqfiles, self.args.outdir)
        else:
            print("Please provide either read1/read2/output files or fastqfiles/outdir")
            sys_exit(1)

    def _run_reads_one(self, read1_filename, read2_filename):
        """Run on a single pair of reads"""
        if read1_filename is not None and read2_filename is not None:
            self._runpairedreads()
        elif read1_filename is not None and read2_filename is None:
            self._runsingleread()
        else:
            self._runnone()

    def _runpairedreads(self):
        """Run on paired reads"""
        self._chk_exists([self.args.read1_filename, self.args.read2_filename])
        top = TopLevel(self.args.config_filename)
        top.run(self.args.output_filename, self.args.read1_filename, self.args.read2_filename)

    def _runsingleread(self):
        """Run on a single read"""
        raise NotImplementedError("Single read processing not implemented")

    def _runnone(self):
        """Handle case with no reads"""
        print("No read files provided")
        sys_exit(1)

    def _run_pairedreads_all(self, fastqfiles, outdir):
        """Run on all paired reads in a directory"""
        if not exists(outdir) or not isdir(outdir):
            print(f'Output directory either does not exist or is not a directory: {outdir}')
            sys_exit(1)
        if not fastqfiles:
            print(f'No fastq files were found with --fastqfiles {fastqfiles}')
            sys_exit(1)

        tic = default_timer()
        readpairs = get_readpairs(fastqfiles)
        nts = get_readpair_files(readpairs)
        len_nts = len(nts)
        
        top = TopLevel(self.args.config_filename)
        for idx, ntd in enumerate(nts, 1):
            hms = timedelta(seconds=default_timer() - tic)
            print(f'hours:mins:secs {hms} -- {idx:4} of {len_nts} readpair {ntd}')
            output = join(outdir, ntd.fcsv)
            try:
                top.run(output, ntd.fread1, ntd.fread2)
            except ValueError as err:
                print(f'Failed analyzing read pair: {ntd.fread1} {ntd.fread2}') 