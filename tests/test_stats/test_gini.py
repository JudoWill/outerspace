"""Tests for Gini coefficient calculation"""

import pytest
from outerspace.stats import GiniCoefficient


def test_gini_perfect_equality(equal_umi):
    """Test Gini coefficient calculation for perfectly equal distribution"""
    result = GiniCoefficient.calculate(equal_umi)
    gini = result['gini_coefficient']
    assert abs(gini) < 0.0001


def test_gini_perfect_inequality(unequal_umi):
    """Test Gini coefficient calculation for perfectly unequal distribution"""
    result = GiniCoefficient.calculate(unequal_umi)
    gini = result['gini_coefficient']
    assert 0.8 < gini < 1.0


def test_gini_moderate_inequality(moderate_umi):
    """Test Gini coefficient calculation for moderately unequal distribution"""
    result = GiniCoefficient.calculate(moderate_umi)
    gini = result['gini_coefficient']
    assert 0 <= gini < 1
    assert 0.2 < gini < 0.3


def test_gini_with_allowed_list(partial_umi):
    """Test Gini coefficient calculation with allowed list including missing keys"""
    # Define allowed list with some missing keys
    allowed_list = ["AAAAAA", "TTTTTT", "CCCCCC", "GGGGGG", "ATATAT", "MISSING1", "MISSING2"]
    
    # Calculate Gini with and without allowed list
    result_without = GiniCoefficient.calculate(partial_umi)
    result_with = GiniCoefficient.calculate(partial_umi, allowed_list=allowed_list)
    gini_without = result_without['gini_coefficient']
    gini_with = result_with['gini_coefficient']
    
    # Gini with allowed list should be higher because it includes zero counts
    assert gini_with > gini_without
    
    # Test with all missing keys
    all_missing = ["MISSING1", "MISSING2", "MISSING3", "MISSING4", "MISSING5"]
    result_all_missing = GiniCoefficient.calculate(partial_umi, allowed_list=all_missing)
    gini_all_missing = result_all_missing['gini_coefficient']
    assert gini_all_missing is None  # No data should return None
    
    # Test with mix of present and missing keys
    mixed_list = ["AAAAAA", "MISSING1", "TTTTTT", "MISSING2", "CCCCCC"]
    result_mixed = GiniCoefficient.calculate(partial_umi, allowed_list=mixed_list)
    gini_mixed = result_mixed['gini_coefficient']
    assert 0 < gini_mixed < 1  # Should be between 0 and 1
    assert gini_mixed > gini_without  # Should be higher than without allowed list


def test_gini_empty_input(empty_umi):
    """Test Gini coefficient calculation with empty input"""
    result = GiniCoefficient.calculate(empty_umi)
    assert result['gini_coefficient'] is None


def test_gini_zero_counts():
    """Test Gini coefficient calculation with all zero counts"""
    from outerspace.umi import UMI
    umi = UMI(mismatches=0)
    sequences = ["AAAAAA", "TTTTTT", "CCCCCC"]
    for seq in sequences:
        umi.consume(seq, 0)  # Add with zero count
    umi.create_mapping()
    
    result = GiniCoefficient.calculate(umi)
    assert result['gini_coefficient'] is None

