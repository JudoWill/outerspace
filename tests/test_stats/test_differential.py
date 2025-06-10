"""Tests for differential abundance calculations"""

import pytest
import numpy as np
import pandas as pd
from outerspace.stats import DifferentialAbundance


def test_differential_abundance_calculation():
    """Test differential abundance calculation with known counts"""
    # Create test DataFrame
    counts_df = pd.DataFrame({
        'sample1': [10, 20, 30, 0],
        'sample2': [20, 10, 0, 30],
        'sample3': [15, 25, 35, 5],
        'sample4': [25, 15, 5, 35],
    }, index=[b"AAAAAA", b"TTTTTT", b"CCCCCC", b"GGGGGG"])
    
    group_map = {
        'sample1': 'group1',
        'sample2': 'group2',
        'sample3': 'group1',
        'sample4': 'group2',
    }
    
    results = DifferentialAbundance.calculate_differential_abundance(counts_df, group_map)
    
    mean_group1 = (10+15)/2
    mean_group2 = (20+25)/2

    fc = mean_group2 / mean_group1
    log2_fc = np.log2(fc)

    #print(fc, log2_fc, mean_group1, mean_group2)
    #print(results.loc[b"AAAAAA", 'log2_fold_change'], results.loc[b"AAAAAA", 'mean_group1'], results.loc[b"AAAAAA", 'mean_group2'])

    # Check results for AAAAAA
    assert abs(results.loc[b"AAAAAA", 'log2_fold_change'] - log2_fc) < 0.01
    assert abs(results.loc[b"AAAAAA", 'mean_group1'] - mean_group1) < 0.01
    assert abs(results.loc[b"AAAAAA", 'mean_group2'] - mean_group2) < 0.01
    assert 0 <= results.loc[b"AAAAAA", 'p_value'] <= 1
    assert results.loc[b"AAAAAA", 'effect_size'] > 0


def test_differential_abundance_identical():
    """Test differential abundance calculation with identical groups"""
    counts_df = pd.DataFrame({
        'sample1': [10, 20, 30],
        'sample2': [10, 20, 30],
        'sample3': [10, 20, 30],
        'sample4': [10, 20, 30],
    }, index=[b"AAAAAA", b"TTTTTT", b"CCCCCC"])
    
    group_map = {
        'sample1': 'group1',
        'sample2': 'group1',
        'sample3': 'group2',
        'sample4': 'group2',
    }
    
    results = DifferentialAbundance.calculate_differential_abundance(counts_df, group_map)
    
    print(results)
    # All log2 fold changes should be 0
    assert all(abs(results['log2_fold_change']) < 0.01)
    # All effect sizes should be 0
    assert all(abs(results['effect_size']) < 0.01)


def test_differential_abundance_empty():
    """Test differential abundance calculation with empty DataFrame"""
    counts_df = pd.DataFrame()
    group_map = {
        'sample1': 'group1',
        'sample2': 'group2',
    }
    
    with pytest.raises(ValueError):
        DifferentialAbundance.calculate_differential_abundance(counts_df, group_map)


def test_differential_abundance_zero_counts():
    """Test differential abundance calculation with zero counts"""
    counts_df = pd.DataFrame({
        'sample1': [0, 0, 0],
        'sample2': [0, 0, 0],
    }, index=[b"AAAAAA", b"TTTTTT", b"CCCCCC"])
    
    group_map = {
        'sample1': 'group1',
        'sample2': 'group2',
    }
    
    results = DifferentialAbundance.calculate_differential_abundance(counts_df, group_map)
    assert all(results['log2_fold_change'] == 0)
    assert all(results['effect_size'] == 0)


def test_differential_abundance_invalid_groups():
    """Test differential abundance calculation with invalid group mapping"""
    counts_df = pd.DataFrame({
        'sample1': [10, 20],
        'sample2': [20, 10],
        'sample3': [15, 25],
    }, index=[b"AAAAAA", b"TTTTTT"])
    
    group_map = {
        'sample1': 'group1',
        'sample2': 'group2',
        'sample3': 'group3',  # Invalid: more than 2 groups
    }
    
    with pytest.raises(ValueError):
        DifferentialAbundance.calculate_differential_abundance(counts_df, group_map)


def test_differential_abundance_with_allowed_list(partial_umi):
    """Test differential abundance calculation with allowed list"""
    allowed_list = ["AAAAAA", "TTTTTT", "CCCCCC", "GGGGGG", "ATATAT"]
    result = DifferentialAbundance.calculate([partial_umi, partial_umi], ['group1', 'group2'], allowed_list=allowed_list)
    assert result is not None
    assert isinstance(result, pd.DataFrame)
    assert all(col in result.columns for col in ['log2_fold_change', 'effect_size', 'p_value']) 