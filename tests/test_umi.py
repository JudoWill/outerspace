"""Tests for UMI clustering and correction"""

import pytest
from outerspace.umi import UMI
from pathlib import Path
import tempfile
import csv
import os


def test_basic_umi_merging():
    """Test basic UMI merging with single mismatch"""
    umi = UMI(mismatches=2)
    
    # Add UMIs that should be merged
    umi.consume("ATCG")
    umi.consume("ATCC")  # One mismatch from ATCG
    umi.consume("GCTA")  # Different UMI
    
    umi.create_mapping()
    
    # Check that similar UMIs are merged
    assert umi["ATCG"] == umi["ATCC"]  # Should be merged
    assert umi["GCTA"] != umi["ATCG"]  # Should not be merged
    
    # Check counts
    counts = umi.corrected_counts
    assert len(counts) == 2  # Should have 2 unique clusters
    assert sum(counts.values()) == 3  # Total count should be 3


def test_pre_clustered_data():
    """Test loading pre-clustered data"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        writer = csv.writer(f)
        writer.writerow(['barcode'])
        writer.writerow(['ATCG'])
        writer.writerow(['ATCC'])
        writer.writerow(['GCTA'])
        temp_file = f.name
    
    try:
        # Load with correction disabled
        umi = UMI.from_csv(temp_file, column='barcode', correct=False)
        
        # Check that UMIs are not merged
        assert umi["ATCG"] != umi["ATCC"]  # Should not be merged
        assert umi["GCTA"] != umi["ATCG"]  # Should not be merged
        
        # Check counts - should be same as original
        counts = umi.corrected_counts
        assert len(counts) == 3  # Should have 3 unique clusters
        assert sum(counts.values()) == 3  # Total count should be 3
        
        # Verify mapping is identity
        assert all(umi[bc] == bc for bc in counts.keys())
        
        # Compare with corrected version
        umi_corrected = UMI.from_csv(temp_file, column='barcode', correct=True)
        assert len(umi_corrected.corrected_counts) < len(counts)  # Should have fewer clusters
        
    finally:
        os.unlink(temp_file)


def test_csv_loading():
    """Test loading UMIs from CSV file"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        writer = csv.writer(f)
        writer.writerow(['barcode'])
        writer.writerow(['ATCG'])
        writer.writerow(['ATCC'])
        writer.writerow(['GCTA'])
        temp_file = f.name
    
    try:
        # Test single file loading
        umi = UMI.from_csv(temp_file, column='barcode', mismatches=2)
        assert umi["ATCG"] == umi["ATCC"]  # Should be merged
        assert umi["GCTA"] != umi["ATCG"]  # Should not be merged
        
        # Test multiple file loading
        temp_dir = tempfile.mkdtemp()
        try:
            # Create multiple CSV files
            for i in range(2):
                with open(os.path.join(temp_dir, f'test_{i}.csv'), 'w') as f:
                    writer = csv.writer(f)
                    writer.writerow(['barcode'])
                    writer.writerow(['ATCG'])
                    writer.writerow(['ATCC'])
            
            umi = UMI.from_csvs(temp_dir, column='barcode', mismatches=1)
            assert umi["ATCG"] == umi["ATCC"]  # Should be merged
            assert umi.corrected_counts[umi["ATCG"]] == 4  # Should have 4 total counts
            
        finally:
            # Cleanup
            for file in Path(temp_dir).glob('*.csv'):
                file.unlink()
            os.rmdir(temp_dir)
            
    finally:
        os.unlink(temp_file)


def test_empty_input():
    """Test handling of empty input"""
    umi = UMI()
    umi.create_mapping()
    assert len(umi.corrected_counts) == 0
    
    # Test empty CSV
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        writer = csv.writer(f)
        writer.writerow(['barcode'])
        temp_file = f.name
    
    try:
        umi = UMI.from_csv(temp_file, column='barcode')
        assert len(umi.corrected_counts) == 0
    finally:
        os.unlink(temp_file)


def test_invalid_input():
    """Test handling of invalid input"""
    umi = UMI()
    
    # Test invalid characters
    with pytest.raises(UnicodeEncodeError):
        umi.consume("ATCG\u2603")  # Snowman character
    
    # Test invalid CSV column
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        writer = csv.writer(f)
        writer.writerow(['wrong_column'])
        writer.writerow(['ATCG'])
        temp_file = f.name
    
    try:
        with pytest.raises(ValueError):
            UMI.from_csv(temp_file, column='barcode')
    finally:
        os.unlink(temp_file)


def test_gini_perfect_equality():
    """Test Gini coefficient calculation for perfectly equal distribution"""
    umi = UMI(mismatches=0)
    
    # Add 10 UMIs with equal counts using sequences that won't merge
    sequences = [
        "AAAAAA", "TTTTTT", "CCCCCC", "GGGGGG",
        "ATATAT", "TATATA", "CGCGCG", "GCGCGC",
        "ACACAC", "TGTGTG"
    ]
    
    for seq in sequences:
        umi.consume(seq)
    
    umi.create_mapping()
    
    # Gini coefficient should be 0 for perfect equality
    assert abs(umi.gini_coefficient()) < 0.0001


def test_gini_perfect_inequality():
    """Test Gini coefficient calculation for perfectly unequal distribution"""
    umi = UMI(mismatches=0)
    
    # Add 10 UMIs with highly unequal distribution
    # One UMI has 100 counts, others have 1
    sequences = [
        "AAAAAA", "TTTTTT", "CCCCCC", "GGGGGG",
        "ATATAT", "TATATA", "CGCGCG", "GCGCGC",
        "ACACAC", "TGTGTG"
    ]
    
    # Add 100 counts of the first sequence
    for _ in range(100):
        umi.consume(sequences[0])
    
    # Add 1 count of each other sequence
    for seq in sequences[1:]:
        umi.consume(seq)
    
    umi.create_mapping()
    
    # Gini coefficient should be close to 1 for perfect inequality
    assert 0.8 < umi.gini_coefficient() < 1.0


def test_gini_moderate_inequality():
    """Test Gini coefficient calculation for moderately unequal distribution"""
    umi = UMI(mismatches=0)
    
    # Add 10 UMIs with moderate inequality
    # Distribution: 1, 2, 3, 4, 5, 5, 4, 3, 2, 1
    sequences = [
        "AAAAAA", "TTTTTT", "CCCCCC", "GGGGGG",
        "ATATAT", "TATATA", "CGCGCG", "GCGCGC",
        "ACACAC", "TGTGTG"
    ]
    counts = [1, 2, 3, 4, 5, 5, 4, 3, 2, 1]
    
    for seq, count in zip(sequences, counts):
        for _ in range(count):
            umi.consume(seq)
    
    umi.create_mapping()
    
    # Gini coefficient should be between 0 and 1
    gini = umi.gini_coefficient()
    assert 0 <= gini < 1
    
    # For this specific distribution, Gini should be around 0.2-0.3
    assert 0.2 < gini < 0.3


def test_gini_with_allowed_list():
    """Test Gini coefficient calculation with allowed list including missing keys"""
    umi = UMI(mismatches=0)
    
    # Add 5 UMIs with unequal distribution
    sequences = ["AAAAAA", "TTTTTT", "CCCCCC", "GGGGGG", "ATATAT"]
    counts = [10, 8, 6, 4, 2]  # Unequal distribution
    
    for seq, count in zip(sequences, counts):
        for _ in range(count):
            umi.consume(seq)
    
    umi.create_mapping()
    
    # Define allowed list with some missing keys
    allowed_list = ["AAAAAA", "TTTTTT", "CCCCCC", "GGGGGG", "ATATAT", "MISSING1", "MISSING2"]
    
    # Calculate Gini with and without allowed list
    gini_without = umi.gini_coefficient()
    gini_with = umi.gini_coefficient(allowed_list=allowed_list)
    print(gini_without, gini_with)
    
    # Gini with allowed list should be higher because it includes zero counts
    assert gini_with > gini_without
    
    # Test with all missing keys
    all_missing = ["MISSING1", "MISSING2", "MISSING3", "MISSING4", "MISSING5"]
    gini_all_missing = umi.gini_coefficient(allowed_list=all_missing)
    assert gini_all_missing is None  # No data should return None
    
    # Test with mix of present and missing keys
    mixed_list = ["AAAAAA", "MISSING1", "TTTTTT", "MISSING2", "CCCCCC"]
    gini_mixed = umi.gini_coefficient(allowed_list=mixed_list)
    assert 0 < gini_mixed < 1  # Should be between 0 and 1
    assert gini_mixed > gini_without  # Should be higher than without allowed list
