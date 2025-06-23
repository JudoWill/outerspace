"""Sequence extraction command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

from argparse import ArgumentParser
from os.path import exists

from outerspace.cli.commands.base import BaseCommand
from outerspace.config import Cfg
from outerspace.pattern import Pattern
from outerspace.read import Read

def run(file_pathname_cfg, file_pathname_fastq1, file_pathname_fastq2=None, filename_csv=None):
    """Parse reads, find user-given patterns, save sequence matches to csv.
    
    This is a convenience wrapper around the findseq command for programmatic use.
    """
    if not exists(file_pathname_cfg):
        raise ValueError(f'Configuration file not found: {file_pathname_cfg}')
    if not exists(file_pathname_fastq1):
        raise ValueError(f'Read 1 file not found: {file_pathname_fastq1}')
    if file_pathname_fastq2 and not exists(file_pathname_fastq2):
        raise ValueError(f'Read 2 file not found: {file_pathname_fastq2}')

    # Create command object
    cmd = FindSeqCommand()
    
    # Set up arguments
    cmd.args = type('Args', (), {
        'config': file_pathname_cfg,
        'read1_filename': file_pathname_fastq1,
        'read2_filename': file_pathname_fastq2,
        'output_filename': filename_csv
    })()

    # Run the command
    cmd.run()

class FindSeqCommand(BaseCommand):
    """Command for extracting sequences from various file formats"""
    
    def _init_parser(self, subparsers):
        """Initialize command-specific argument parser"""
        parser = subparsers.add_parser('findseq',
            help='Extract sequences from files based on configuration patterns')
        parser.add_argument('config',
            help='Configuration file with search patterns')
        parser.add_argument('-1', '--read1_filename',
            help='Input file for read 1 (FASTQ, FASTA, SAM, BAM) or single read file')
        parser.add_argument('-2', '--read2_filename',
            help='Input file for read 2 (FASTQ, FASTA, SAM, BAM) for paired reads')
        parser.add_argument('-o', '--output_filename',
            help='Output CSV file name')
        parser.add_argument('--region',
            help='SAM/BAM region specification (e.g., "chr1:1-1000")')
        parser.add_argument('--fetch',
            help='SAM/BAM fetch mode (mapped, unmapped, all)')
        return parser

    def _detect_file_format(self, filename):
        """Auto-detect file format based on extension"""
        key = filename.lower().replace('.gz', '')
        if key.endswith(('.fastq', '.fq')):
            return 'fastq'
        elif key.endswith(('.fasta', '.fa')):
            return 'fasta'
        elif key.endswith('.sam'):
            return 'sam'
        elif key.endswith('.bam'):
            return 'bam'
        else:
            raise ValueError(f"Cannot detect format for file: {filename}. Please use a file with .fastq, .fasta, .sam, or .bam extension.")

    def _process_single_file(self, filename):
        """Process a single file of any supported format"""
        format_type = self._detect_file_format(filename)
        
        if format_type in ['sam', 'bam']:
            # Handle SAM/BAM with optional region and fetch parameters
            fetch_param = None
            if self.args.region:
                fetch_param = self.args.region
            elif self.args.fetch:
                if self.args.fetch == 'mapped':
                    fetch_param = True
                elif self.args.fetch == 'unmapped':
                    fetch_param = False
                elif self.args.fetch == 'all':
                    fetch_param = None
                else:
                    raise ValueError(f"Invalid fetch mode: {self.args.fetch}")
            
            return Read.from_bam(filename, fetch=fetch_param)
        elif format_type == 'fasta':
            return Read.from_fasta(filename)
        elif format_type == 'fastq':
            return Read.from_fastq(filename)
        else:
            raise ValueError(f"Unsupported format: {format_type}")

    def _process_paired_files(self, file1, file2):
        """Process paired files of any supported format"""
        format1 = self._detect_file_format(file1)
        format2 = self._detect_file_format(file2)
        
        if format1 != format2:
            raise ValueError(f"Paired files must have the same format. Got {format1} and {format2}")
        
        return Read.from_paired_fastx(file1, file2, format=format1)

    def _process_reads_with_patterns(self, reads):
        """Process reads using Pattern objects and return results"""
        patterns = self._parse_patterns_from_config()
        results = []
        
        for read in reads:
            read_results = {}
            for pattern in patterns:
                hits = pattern.search(read)
                if hits:
                    # Handle multiple, first, last logic
                    if pattern.multiple == 'first':
                        hits = [hits] if not isinstance(hits, list) else [hits[0]]
                    elif pattern.multiple == 'last':
                        hits = [hits] if not isinstance(hits, list) else [hits[-1]]
                    # hits is already a list for 'all'
                    
                    for hit in hits:
                        read_results.update(hit.captured)
            
            if read_results:
                read_results['read_id'] = read.name
                results.append(read_results)
        
        return results

    def _write_results_to_csv(self, results, output_filename):
        """Write results to CSV file"""
        import csv
        
        if not results:
            print("No matches found")
            return
        
        # Get all unique keys from all results
        all_keys = set()
        for result in results:
            all_keys.update(result.keys())
        
        # Sort keys for consistent output
        fieldnames = sorted(all_keys)
        
        with open(output_filename, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
        
        print(f"Results written to {output_filename}")

    def run(self):
        """Run the findseq command"""
        # Load config if provided
        if self.args.config:
            self._load_config(self.args.config)
        
        # Merge config and args with defaults
        defaults = {
            'read1_filename': None,
            'read2_filename': None,
            'output_filename': None,
            'region': None,
            'fetch': None
        }
        self._merge_config_and_args(defaults)

        if self.args.config is None:
            raise ValueError('Please provide a config filename')

        # Check for required output file
        if not self.args.output_filename:
            raise ValueError('Please provide an output filename with -o')

        # Process based on input files
        if self.args.read1_filename and self.args.read2_filename:
            # Paired files
            self._chk_exists([self.args.read1_filename, self.args.read2_filename])
            reads = self._process_paired_files(self.args.read1_filename, self.args.read2_filename)
        elif self.args.read1_filename:
            # Single file
            self._chk_exists([self.args.read1_filename])
            reads = self._process_single_file(self.args.read1_filename)
        else:
            raise ValueError("Please provide either -1 for single file or -1/-2 for paired files")

        # Process reads with patterns
        results = self._process_reads_with_patterns(reads)
        
        # Write results
        self._write_results_to_csv(results, self.args.output_filename) 