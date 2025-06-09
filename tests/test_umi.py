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
