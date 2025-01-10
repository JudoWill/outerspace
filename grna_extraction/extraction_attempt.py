"""Extraction of sequence from reads"""
####
####__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
####__author__ = "SC Barrera"
####
####from os.path import join
####from Bio import SeqIO
####import gzip
####
####from itertools import islice
####from Bio.Seq import reverse_complement
####
####import regex
####import tqdm
####import csv
####
####dir_crispr = '../../../nonn-lab/rachel-test-crispr'
####
####dir_read = join(dir_crispr, 'reads/')
####path1 = join(dir_read, '409-4_S1_L001_R1_001.fastq.gz')
####path2 = join(dir_read, '409-4_S1_L001_R2_001.fastq.gz')
####
####def iterate_reads(path):
####    "Iterate reads from a gzipped or regular fastq files"
####    
####    if path.endswith('.gz'):
####        with gzip.open(path, mode='rt') as handle:
####            for read in SeqIO.parse(handle, 'fastq'):
####                yield read
####    else:
####        with open(path) as handle:
####            for read in SeqIO.parse(handle, 'fastq'):
####                yield read
####
####def iterate_readpairs(path1, path2):
####    
####    reads1 = iterate_reads(path1)
####    reads2 = iterate_reads(path2)
####    
####    for r1, r2 in zip(reads1, reads2):
####        yield r1, r2
####                
####def get_capture_from_read(reg_exp, read):
####    
####    sequence = str(read.seq)
####    result = reg_exp.findall(sequence)
####    
####    if len(result) == 1:
####        return result[0]
####    
####def extract_from_paired_reads(read1_path, read2_path, forward_reg, proto_reg, reverse_reg):
####    
####    total = 0
####    missed = 0
####    
####    for read1, read2 in tqdm.tqdm(iterate_readpairs(path1, path2)):
####        total += 1
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
####    print(f'Total sequences: {total}\nMissed sequences: {missed}')
####    
def process_paired_read_file(outpath, path1, path2, forward_umi_reg, protospacer_reg, reverse_umi_reg):
    print('YAYY')
####    
####    with open(outpath, mode='w') as handle:
####        fieldnames = ['read_id', 'forward_umi', 'protospacer', 'reverse_umi']
####        writer = csv.DictWriter(handle, fieldnames)
####        writer.writeheader()
####        
####        stream = extract_from_paired_reads(path1, path2, forward_umi_reg, protospacer_reg, reverse_umi_reg)
####        writer.writerows(stream)
####
####forward_umi_reg = regex.compile('(?P<UMI>.{8})(?:CTTGGCTTTATATATCTTGTGG){s<=4}', flags=regex.BESTMATCH)
####protospacer_reg = regex.compile('(?:TATCTTGTGGAAAGGACGAAACACC){s<=4}(?P<protospacer>.{19,21})(?:GTTTAAGTACTCTGTGCTGGAAACAG){s<=4}',
####                                flags=regex.BESTMATCH)
####
####back_umi_forward = 'gtgtgtcagttagggtgtggaa'.upper()
####back_umi_rc = reverse_complement(back_umi_forward)
####
####reverse_umi_reg = regex.compile(f'(?P<UMI>.{{8}})(?:{back_umi_rc}){{s<=4}}', flags=regex.BESTMATCH)
####
##### Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
