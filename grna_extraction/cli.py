"""Extraction of sequence from reads"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "RE Berman"

from sys import exit as sys_exit
from os.path import isdir
from os.path import isfile
from os.path import exists
from argparse import ArgumentParser
from logging import debug
from grna_extraction.config import Cfg

class Cli:
    """The command line interface argument requirements"""
    def __init__(self, args=None):
        self.args = self._init_args(args)
        self.cfg = self._init_config()
        
    @staticmethod
    def _init_args(args=None):
        """cli for entering the reads"""
        parser = ArgumentParser(
                        prog='grna_extraction',
                        description='Get protospacers and UMIs',
                        epilog='Created by ThreeBlindMice - See how they run code')
        parser.add_argument('config_filename', nargs='?',
            help='Configuration file holds user input search patterns ')
        parser.add_argument('read1_filename', nargs='*',
            help='zipped fastq files or a directory containing fastq files')
        parser.add_argument('read2_filename', nargs='*',
            help='zipped fastq files or a directory containing fastq files')
        parser.add_argument('output_filename', nargs='*',
            help='zipped fastq files or a directory containing fastq files')

        args = parser.parse_args(args)
        print(f'ARGS: {args}')
        debug(f'ARGS: {args}') 
        #if args.config_filename:
            #_chk_exists(args.config_filename)
        #else:
            #print('EXITING: NEED FASTQ READS FILES TO READ')
            #sys_exit(1)
        return args

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

    def _init_config(self):
        """checking that config exists"""
        filenamecfg = self.args.config_filename
        if filenamecfg is None:
            # TODO: give researcher explicit instructions how to generate a config file - give hints or examples
            print("Provide configuration file name")
            return None
        if exists(filenamecfg):
            return Cfg(filenamecfg)
        print(f'Configuration file does not exist: {filenamecfg}')
        return None

# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
