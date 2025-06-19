"""Tests for Read class"""

import pytest
import tempfile
import os
from outerspace.read import Read
import pysam


def test_read_creation():
    """Test basic Read object creation"""
    read = Read(seq="ATCG", pair="R1", name="test_read")
    
    assert read.seq == "ATCG"
    assert read.pair == "R1"
    assert read.name == "test_read"


def test_read_creation_without_name():
    """Test Read object creation without name parameter"""
    read = Read(seq="ATCG", pair="R1")
    
    assert read.seq == "ATCG"
    assert read.pair == "R1"
    assert read.name is None


def test_reverse_complement():
    """Test reverse complement property"""
    read = Read(seq="ATCG", pair="R1")
    
    # Test basic reverse complement
    assert read.seq_rc == "CGAT"
    
    # Test with different sequence
    read2 = Read(seq="GCTA", pair="R1")
    assert read2.seq_rc == "TAGC"
    
    # Test with longer sequence
    read3 = Read(seq="ATCGATCG", pair="R1")
    assert read3.seq_rc == "CGATCGAT"


def test_reverse_complement_edge_cases():
    """Test reverse complement with edge cases"""
    # Test empty sequence
    read = Read(seq="", pair="R1")
    assert read.seq_rc == ""
    
    # Test single nucleotide
    read = Read(seq="A", pair="R1")
    assert read.seq_rc == "T"
    
    # Test with N (should remain N)
    read = Read(seq="ATCN", pair="R1")
    assert read.seq_rc == "NGAT"


def test_reverse_complement_round_trip():
    """Test that reverse complement of reverse complement equals original"""
    sequences = ["ATCG", "GCTA", "ATCGATCG", "N", ""]
    
    for seq in sequences:
        read = Read(seq=seq, pair="R1")
        # RC of RC should equal original
        assert read.seq_rc[::-1].translate(str.maketrans('ACGT', 'TGCA')) == seq


def test_string_representation():
    """Test string and repr methods"""
    read = Read(seq="ATCG", pair="R1", name="test")
    
    # Test __str__
    assert str(read) == "ATCG\nR1"
    
    # Test __repr__
    assert repr(read) == "ATCG\nR1"


def test_string_representation_without_name():
    """Test string representation when name is None"""
    read = Read(seq="ATCG", pair="R1")
    
    assert str(read) == "ATCG\nR1"
    assert repr(read) == "ATCG\nR1"


def test_sequence_properties():
    """Test various sequence properties"""
    read = Read(seq="ATCG", pair="R1")
    
    # Test sequence length
    assert len(read.seq) == 4
    
    # Test sequence content
    assert "A" in read.seq
    assert "T" in read.seq
    assert "C" in read.seq
    assert "G" in read.seq


def test_pair_assignment():
    """Test different pair assignments"""
    # Test R1
    read1 = Read(seq="ATCG", pair="R1")
    assert read1.pair == "R1"
    
    # Test R2
    read2 = Read(seq="ATCG", pair="R2")
    assert read2.pair == "R2"
    
    # Test other pair names
    read3 = Read(seq="ATCG", pair="both")
    assert read3.pair == "both"


def test_sequence_mutation():
    """Test that sequence can be modified"""
    read = Read(seq="ATCG", pair="R1")
    
    # Modify sequence
    read.seq = "GCTA"
    assert read.seq == "GCTA"
    assert read.seq_rc == "TAGC"


def test_name_assignment():
    """Test name assignment and modification"""
    read = Read(seq="ATCG", pair="R1")
    
    # Initially None
    assert read.name is None
    
    # Assign name
    read.name = "new_name"
    assert read.name == "new_name"
    
    # Create with name
    read2 = Read(seq="ATCG", pair="R1", name="initial_name")
    assert read2.name == "initial_name"


def test_sequence_validation():
    """Test sequence validation (basic DNA characters)"""
    # Valid sequences
    valid_seqs = ["ATCG", "GCTA", "ATCGATCG", "N", "ATCN"]
    
    for seq in valid_seqs:
        read = Read(seq=seq, pair="R1")
        assert read.seq == seq
    
    # Test with lowercase (should work)
    read = Read(seq="atcg", pair="R1")
    assert read.seq == "atcg"
    assert read.seq_rc == "cgat"


def test_empty_sequence():
    """Test handling of empty sequences"""
    read = Read(seq="", pair="R1")
    
    assert read.seq == ""
    assert read.seq_rc == ""
    assert len(read.seq) == 0


def test_sequence_with_spaces():
    """Test handling of sequences with spaces"""
    read = Read(seq=" ATCG ", pair="R1")
    
    assert read.seq == " ATCG "
    # Reverse complement should preserve spaces
    assert read.seq_rc == " CGAT "


def test_multiple_reads():
    """Test creating multiple read objects"""
    reads = [
        Read(seq="ATCG", pair="R1", name="read1"),
        Read(seq="GCTA", pair="R2", name="read2"),
        Read(seq="ATCGATCG", pair="R1", name="read3")
    ]
    
    assert len(reads) == 3
    assert reads[0].seq == "ATCG"
    assert reads[1].seq == "GCTA"
    assert reads[2].seq == "ATCGATCG"
    
    assert reads[0].pair == "R1"
    assert reads[1].pair == "R2"
    assert reads[2].pair == "R1"


def test_reverse_complement_consistency():
    """Test that reverse complement is consistent across multiple calls"""
    read = Read(seq="ATCG", pair="R1")
    
    # Multiple calls should return same result
    rc1 = read.seq_rc
    rc2 = read.seq_rc
    rc3 = read.seq_rc
    
    assert rc1 == rc2 == rc3 == "CGAT"


def test_sequence_immutability():
    """Test that sequence property behaves correctly"""
    read = Read(seq="ATCG", pair="R1")
    
    # Get sequence
    original_seq = read.seq
    original_rc = read.seq_rc
    
    # Modify sequence
    read.seq = "GCTA"
    
    # Original variables should not change
    assert original_seq == "ATCG"
    assert original_rc == "CGAT"
    
    # But read properties should change
    assert read.seq == "GCTA"
    assert read.seq_rc == "TAGC"


# File I/O Tests

def test_from_fasta():
    """Test reading from FASTA file"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(">read1\nATCG\n")
        f.write(">read2\nGCTA\n")
        f.write(">read3\nATCGATCG\n")
        temp_file = f.name
    
    try:
        reads = list(Read.from_fasta(temp_file, read_pair='R1'))
        
        assert len(reads) == 3
        assert reads[0].seq == "ATCG"
        assert reads[0].name == "read1"
        assert reads[0].pair == "R1"
        
        assert reads[1].seq == "GCTA"
        assert reads[1].name == "read2"
        assert reads[1].pair == "R1"
        
        assert reads[2].seq == "ATCGATCG"
        assert reads[2].name == "read3"
        assert reads[2].pair == "R1"
        
    finally:
        os.unlink(temp_file)


def test_from_fastq():
    """Test reading from FASTQ file"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
        f.write("@read1\nATCG\n+\nIIII\n")
        f.write("@read2\nGCTA\n+\nIIII\n")
        f.write("@read3\nATCGATCG\n+\nIIIIIIII\n")
        temp_file = f.name
    
    try:
        reads = list(Read.from_fastq(temp_file, read_pair='R2'))
        
        assert len(reads) == 3
        assert reads[0].seq == "ATCG"
        assert reads[0].name == "read1"
        assert reads[0].pair == "R2"
        
        assert reads[1].seq == "GCTA"
        assert reads[1].name == "read2"
        assert reads[1].pair == "R2"
        
        assert reads[2].seq == "ATCGATCG"
        assert reads[2].name == "read3"
        assert reads[2].pair == "R2"
        
    finally:
        os.unlink(temp_file)


def test_from_paired_fastx_fasta():
    """Test reading from paired FASTA files"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='_1.fasta', delete=False) as f1:
        f1.write(">read1\nATCG\n")
        f1.write(">read2\nGCTA\n")
        temp_file1 = f1.name
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='_2.fasta', delete=False) as f2:
        f2.write(">read1\nTAGC\n")
        f2.write(">read2\nCGAT\n")
        temp_file2 = f2.name
    
    try:
        reads = list(Read.from_paired_fastx(temp_file1, temp_file2, format='fasta'))
        
        assert len(reads) == 4  # 2 reads from each file
        
        # Check R1 reads
        assert reads[0].seq == "ATCG"
        assert reads[0].name == "read1"
        assert reads[0].pair == "R1"
        
        assert reads[2].seq == "GCTA"
        assert reads[2].name == "read2"
        assert reads[2].pair == "R1"
        
        # Check R2 reads
        assert reads[1].seq == "TAGC"
        assert reads[1].name == "read1"
        assert reads[1].pair == "R2"
        
        assert reads[3].seq == "CGAT"
        assert reads[3].name == "read2"
        assert reads[3].pair == "R2"
        
    finally:
        os.unlink(temp_file1)
        os.unlink(temp_file2)


def test_from_paired_fastx_fastq():
    """Test reading from paired FASTQ files"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='_1.fastq', delete=False) as f1:
        f1.write("@read1\nATCG\n+\nIIII\n")
        f1.write("@read2\nGCTA\n+\nIIII\n")
        temp_file1 = f1.name
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='_2.fastq', delete=False) as f2:
        f2.write("@read1\nTAGC\n+\nIIII\n")
        f2.write("@read2\nCGAT\n+\nIIII\n")
        temp_file2 = f2.name
    
    try:
        reads = list(Read.from_paired_fastx(temp_file1, temp_file2, format='fastq'))
        
        assert len(reads) == 4  # 2 reads from each file
        
        # Check R1 reads
        assert reads[0].seq == "ATCG"
        assert reads[0].name == "read1"
        assert reads[0].pair == "R1"
        
        assert reads[2].seq == "GCTA"
        assert reads[2].name == "read2"
        assert reads[2].pair == "R1"
        
        # Check R2 reads
        assert reads[1].seq == "TAGC"
        assert reads[1].name == "read1"
        assert reads[1].pair == "R2"
        
        assert reads[3].seq == "CGAT"
        assert reads[3].name == "read2"
        assert reads[3].pair == "R2"
        
    finally:
        os.unlink(temp_file1)
        os.unlink(temp_file2)


def test_from_paired_fastx_invalid_format():
    """Test error handling for invalid format in paired FASTX"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='_1.fasta', delete=False) as f1:
        f1.write(">read1\nATCG\n")
        temp_file1 = f1.name
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='_2.fasta', delete=False) as f2:
        f2.write(">read1\nTAGC\n")
        temp_file2 = f2.name
    
    try:
        with pytest.raises(ValueError, match="Unsupported format"):
            list(Read.from_paired_fastx(temp_file1, temp_file2, format='invalid'))
    finally:
        os.unlink(temp_file1)
        os.unlink(temp_file2)


def test_from_sam():
    """Test reading from SAM file"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.sam', delete=False) as f:
        # Write SAM header
        f.write("@HD\tVN:1.0\tSO:unsorted\n")
        f.write("@SQ\tSN:chr1\tLN:1000\n")
        
        # Write SAM records
        # read1: first in pair (flag 64), mapped
        f.write("read1\t64\tchr1\t1\t255\t4M\t*\t0\t0\tATCG\tIIII\n")
        # read2: second in pair (flag 128), mapped  
        f.write("read2\t128\tchr1\t1\t255\t4M\t*\t0\t0\tGCTA\tIIII\n")
        # read3: single end (no pair flags), mapped
        f.write("read3\t0\tchr1\t1\t255\t8M\t*\t0\t0\tATCGATCG\tIIIIIIII\n")
        # read4: unmapped (flag 4)
        f.write("read4\t4\t*\t0\t0\t*\t*\t0\t0\tATCG\tIIII\n")
        temp_file = f.name
    
    try:
        # Test without fetch parameter - should get all reads
        reads = list(Read.from_bam(temp_file))
        
        assert len(reads) == 4  # All reads including unmapped
        
        # Check read1 (first in pair)
        assert reads[0].seq == "ATCG"
        assert reads[0].name == "read1"
        assert reads[0].pair == "R1"
        
        # Check read2 (second in pair)
        assert reads[1].seq == "GCTA"
        assert reads[1].name == "read2"
        assert reads[1].pair == "R2"
        
        # Check read3 (single end)
        assert reads[2].seq == "ATCGATCG"
        assert reads[2].name == "read3"
        assert reads[2].pair == "R1"  # Default for single end
        
        # Check unmapped read
        assert reads[3].seq == "ATCG"
        assert reads[3].name == "read4"
        assert reads[3].pair == "R1"  # Default for unmapped
        
    finally:
        os.unlink(temp_file)


def test_from_bam():
    """Test reading from BAM file (requires pysam to write BAM)"""
    # Create a simple BAM file using pysam
    with tempfile.NamedTemporaryFile(suffix='.bam', delete=False) as f:
        temp_bam = f.name
    
    try:
        # Create BAM file with pysam
        header = {'HD': {'VN': '1.0'}, 'SQ': [{'LN': 1000, 'SN': 'chr1'}]}
        with pysam.AlignmentFile(temp_bam, "wb", header=header) as bam:
            # Create read1 (first in pair)
            read1 = pysam.AlignedSegment()
            read1.query_name = "read1"
            read1.query_sequence = "ATCG"
            read1.flag = 64  # First in pair
            read1.reference_id = 0
            read1.reference_start = 0
            read1.mapping_quality = 255
            read1.cigarstring = "4M"
            read1.query_qualities = pysam.qualitystring_to_array("IIII")
            bam.write(read1)
            
            # Create read2 (second in pair)
            read2 = pysam.AlignedSegment()
            read2.query_name = "read2"
            read2.query_sequence = "GCTA"
            read2.flag = 128  # Second in pair
            read2.reference_id = 0
            read2.reference_start = 0
            read2.mapping_quality = 255
            read2.cigarstring = "4M"
            read2.query_qualities = pysam.qualitystring_to_array("IIII")
            bam.write(read2)
        
        # Read the BAM file
        reads = list(Read.from_bam(temp_bam))
        
        assert len(reads) == 2
        
        # Check read1 (first in pair)
        assert reads[0].seq == "ATCG"
        assert reads[0].name == "read1"
        assert reads[0].pair == "R1"
        
        # Check read2 (second in pair)
        assert reads[1].seq == "GCTA"
        assert reads[1].name == "read2"
        assert reads[1].pair == "R2"
        
    finally:
        os.unlink(temp_bam)


def test_from_bam_with_fetch_region():
    """Test reading from BAM file with specific region fetch"""
    # Create a simple BAM file using pysam
    with tempfile.NamedTemporaryFile(suffix='.bam', delete=False) as f:
        temp_bam = f.name
    
    try:
        # Create BAM file with pysam
        header = {'HD': {'VN': '1.0'}, 'SQ': [{'LN': 1000, 'SN': 'chr1'}]}
        with pysam.AlignmentFile(temp_bam, "wb", header=header) as bam:
            # Create read1 in region chr1:1-10
            read1 = pysam.AlignedSegment()
            read1.query_name = "read1"
            read1.query_sequence = "ATCG"
            read1.flag = 64
            read1.reference_id = 0
            read1.reference_start = 1
            read1.mapping_quality = 255
            read1.cigarstring = "4M"
            read1.query_qualities = pysam.qualitystring_to_array("IIII")
            bam.write(read1)
            
            # Create read2 in region chr1:50-60
            read2 = pysam.AlignedSegment()
            read2.query_name = "read2"
            read2.query_sequence = "GCTA"
            read2.flag = 128
            read2.reference_id = 0
            read2.reference_start = 50
            read2.mapping_quality = 255
            read2.cigarstring = "4M"
            read2.query_qualities = pysam.qualitystring_to_array("IIII")
            bam.write(read2)
        
        # Index the BAM file for region fetching
        pysam.index(temp_bam)
        
        # Read only reads in region chr1:1-10
        reads = list(Read.from_bam(temp_bam, fetch="chr1:1-10"))
        
        assert len(reads) == 1
        assert reads[0].seq == "ATCG"
        assert reads[0].name == "read1"
        assert reads[0].pair == "R1"
        
        # Test with different region format
        reads = list(Read.from_bam(temp_bam, fetch="chr1:50-60"))
        
        assert len(reads) == 1
        assert reads[0].seq == "GCTA"
        assert reads[0].name == "read2"
        assert reads[0].pair == "R2"
        
        # Test with region that includes both reads
        reads = list(Read.from_bam(temp_bam, fetch="chr1:1-100"))
        
        assert len(reads) == 2
        assert reads[0].seq == "ATCG"
        assert reads[0].name == "read1"
        assert reads[0].pair == "R1"
        assert reads[1].seq == "GCTA"
        assert reads[1].name == "read2"
        assert reads[1].pair == "R2"
        
    finally:
        # Clean up BAM file and index
        if os.path.exists(temp_bam):
            os.unlink(temp_bam)
        index_file = temp_bam + ".bai"
        if os.path.exists(index_file):
            os.unlink(index_file) 