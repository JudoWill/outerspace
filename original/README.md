# GRNA Extraction and Analysis Tools

A collection of Python scripts for processing and analyzing gRNA barcode data, with a focus on UMI (Unique Molecular Identifier) handling and statistical analysis.

## Project Evolution

The project has evolved through several key stages:

1. **Initial Setup** (a1024da, March 31, 2025 08:21:08 -0400)
   - Added `regex` to requirements
   - Set up virtual environment and test targets

2. **Core Functionality** (10dff43 - 70e661c, March 31, 2025 11:03:22 -0400 to March 31, 2025 11:03:49 -0400)
   - Added script to extract and normalize barcodes
   - Implemented counting of barcodes per protospacer
   - Added support for allowed lists in counting
   - Modified scripts to handle directory inputs/outputs

3. **Advanced Features** (c42a621 - 559703a, March 31, 2025 13:39:31 -0400 to March 31, 2025 17:20:08 -0400)
   - Added downsampling capability
   - Implemented Gini coefficient calculation
   - Refactored UMI logic into a dedicated class
   - Updated count.py and collapse.py to use the new UMI class

4. **Metrics and Analysis** (b73c7be - 3acf198, April 1, 2025 09:00:09 -0400 to April 1, 2025 12:19:24 -0400)
   - Added comprehensive metrics
   - Implemented Gini coefficient tests
   - Added barcode and key-column Gini calculations
   - Created a dedicated tool for Gini value calculation

## Scripts



## UMI Class

The `UMI` class (in `grna_extraction/umi.py`) provides core functionality for UMI handling:

### Key Features
- UMI clustering with configurable mismatch tolerance
- Multiple clustering methods (cluster, adjacency, directional)
- Gini coefficient calculation
- CSV file integration
- Pre-clustered data support

### Main Methods
- `consume(umi)`: Add a UMI to the counts
- `create_mapping()`: Create mapping between original and corrected barcodes
- `gini_coefficient()`: Calculate Gini coefficient for UMI distribution
- `from_csv()`: Create UMI object from CSV file
- `from_csvs()`: Create UMI object from multiple CSV files

### Properties
- `corrected_counts`: Get counts of corrected barcodes
- `_counts`: Original UMI counts
- `_mapping`: Mapping between original and corrected barcodes

## Requirements

- Python 3.x
- umi_tools
- matplotlib
- numpy
- tqdm
- yaml

## Usage

Each script includes detailed help documentation accessible via the `--help` flag. For example:

```bash
python collapse.py --help
python count.py --help
python gini.py --help
python visualize.py --help
```

## Testing

The project includes a comprehensive test suite in the `tests` directory, with particular focus on UMI functionality and Gini coefficient calculations. The main test file `test_umi.py` includes the following tests:

### UMI Clustering Tests

```python
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
```

This test verifies the core UMI clustering functionality. It's crucial because it ensures that:
1. Similar sequences (within mismatch tolerance) are correctly merged
2. Dissimilar sequences remain separate
3. Counts are properly aggregated after merging
4. The clustering algorithm correctly handles the specified mismatch threshold

```python
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
```

This test ensures proper handling of pre-clustered data, which is important because:
1. Some datasets may already have undergone UMI correction
2. Users need the ability to skip clustering when appropriate
3. The system must maintain data integrity when clustering is disabled
4. It verifies that the `correct` parameter works as expected

### Input Validation Tests

```python
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
```

This test verifies proper handling of edge cases, which is critical because:
1. Empty inputs are common in real-world scenarios
2. The system must handle them gracefully without errors
3. It ensures no false positives in UMI detection
4. It validates the robustness of the CSV loading functionality

```python
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
```

This test ensures proper error handling, which is important because:
1. Invalid inputs must be caught and handled appropriately
2. Users need clear error messages for debugging
3. The system must maintain data integrity when errors occur
4. It verifies that input validation works as expected

### Gini Coefficient Tests

```python
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
```

This test verifies the Gini coefficient calculation for perfect equality, which is important because:
1. It serves as a baseline for the Gini calculation
2. It ensures the coefficient correctly identifies uniform distributions
3. It validates the mathematical correctness of the implementation
4. It provides a reference point for other inequality tests

```python
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
```

This test verifies the Gini coefficient calculation for perfect inequality, which is important because:
1. It tests the other extreme of the Gini coefficient range
2. It ensures the coefficient correctly identifies highly skewed distributions
3. It validates the sensitivity of the calculation to extreme cases
4. It provides a reference point for moderate inequality tests

```python
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
    
    # Gini with allowed list should be higher because it includes zero counts
    assert gini_with > gini_without
    
    # Test with all missing keys
    all_missing = ["MISSING1", "MISSING2", "MISSING3", "MISSING4", "MISSING5"]
    gini_all_missing = umi.gini_coefficient(allowed_list=all_missing)
    assert gini_all_missing is None  # No data should return None
```

This test verifies the Gini coefficient calculation with allowed lists, which is important because:
1. It tests the handling of missing values in the distribution
2. It ensures the coefficient correctly accounts for zero counts
3. It validates the behavior with user-specified allowed lists
4. It verifies proper handling of edge cases with all missing keys

These tests ensure the reliability of the UMI clustering and Gini coefficient calculations, which are critical for accurate barcode analysis and distribution measurements.

