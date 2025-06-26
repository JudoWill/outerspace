"""Hold read information"""

import pyfastx
import pysam
from itertools import zip_longest

class Read:
    """Hold read information"""

    def __init__(self, seq: str, pair: str, name: str = None):
        self.name = name
        self.seq = seq
        self.pair = pair

    @property
    def seq_rc(self):
        return self.seq[::-1].translate(str.maketrans('ACGTacgt', 'TGCAtgca'))
    
    def __str__(self):
        return f'{self.seq}\n{self.pair}'
    
    def __repr__(self):
        return f'{self.seq}\n{self.pair}'
    
    @classmethod
    def from_fasta(cls, fasta_file: str, read_pair: str = 'R1'):
        """Iterate Read objects from a FASTA file"""
        for name, seq in pyfastx.Fasta(fasta_file, build_index=False):
            yield cls(seq=str(seq), name=name, pair=read_pair)
    
    @classmethod
    def from_fastq(cls, fastq_file: str, read_pair: str = 'R1'):
        """Iterate Read objects from a FASTQ file"""    
        for name, seq, _ in pyfastx.Fastq(fastq_file, build_index=False):
            yield cls(seq=str(seq), name=name, pair=read_pair)
    
    @classmethod
    def from_paired_fastx(cls, fastx_file1: str, fastx_file2: str, format: str = 'fasta'):
        """Iterate Read objects from paired FASTX files"""
        if format.lower() == 'fasta':
            fx1 = Read.from_fasta(fastx_file1, read_pair='R1')
            fx2 = Read.from_fasta(fastx_file2, read_pair='R2')
        elif format.lower() == 'fastq':
            fx1 = Read.from_fastq(fastx_file1, read_pair='R1')
            fx2 = Read.from_fastq(fastx_file2, read_pair='R2')
        else:
            raise ValueError(f"Unsupported format: {format}")
        
        # Iterate through both files in alternating fashion
        for read1, read2 in zip_longest(fx1, fx2):
            if read1 is not None: yield read1
            if read2 is not None: yield read2

    
    @classmethod
    def from_bam(cls, bam_file: str, fetch: str = None):
        """Iterate Read objects from a SAM/BAM file"""
        # Determine file format and open accordingly
        if bam_file.endswith('.sam'):
            mode = "r"  # SAM file
        else:
            mode = "rb"  # BAM file
            
        with pysam.AlignmentFile(bam_file, mode) as bam:
            if fetch is not None:
                # Parse region string like "chr1:1-10" or "chr1"
                if ':' in fetch:
                    # Region with coordinates: "chr1:1-10"
                    contig, coords = fetch.split(':', 1)
                    if '-' in coords:
                        start_str, stop_str = coords.split('-', 1)
                        start = int(start_str) - 1  # Convert to 0-based
                        stop = int(stop_str)  # Keep as 1-based for pysam
                    else:
                        start = int(coords) - 1  # Single position
                        stop = start + 1
                else:
                    # Just contig name: "chr1"
                    contig = fetch
                    start = None
                    stop = None
                
                reads = bam.fetch(contig, start, stop)
            else:
                reads = bam

            for read in reads:
                # Infer read pair from SAM flag
                # Bit 6 (64) indicates first in pair, bit 7 (128) indicates second in pair
                if read.flag & 64:  # First in pair
                    read_pair = 'R1'
                elif read.flag & 128:  # Second in pair
                    read_pair = 'R2'
                else:  # Single end or unpaired
                    read_pair = 'R1'
                
                yield cls(seq=read.query_sequence, name=read.query_name, pair=read_pair)
            