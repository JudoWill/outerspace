"""Tests for single sample differential abundance calculations"""

import pytest
import numpy as np
import pandas as pd
from outerspace.stats.differential import SingleSampleDifferentialAbundance


def test_single_sample_differential_abundance():
    """Test single sample differential abundance calculation with known counts"""
    # Create test DataFrames for each group
    group1_df = pd.DataFrame(
        {
            "sample1": [10, 20, 30, 0],
        },
        index=[b"AAAAAA", b"TTTTTT", b"CCCCCC", b"GGGGGG"],
    )

    group2_df = pd.DataFrame(
        {
            "sample2": [20, 10, 0, 30],
        },
        index=[b"AAAAAA", b"TTTTTT", b"CCCCCC", b"GGGGGG"],
    )

    results = SingleSampleDifferentialAbundance.calculate_differential_abundance_single(
        group1_df, group2_df
    )

    # Check results for AAAAAA
    assert abs(results.loc[b"AAAAAA", "log2_fold_change"] - np.log2(20 / 10)) < 0.01
    assert abs(results.loc[b"AAAAAA", "effect_size"] - 10) < 0.01  # Absolute difference
    assert abs(results.loc[b"AAAAAA", "mean_group1"] - 10) < 0.01
    assert abs(results.loc[b"AAAAAA", "mean_group2"] - 20) < 0.01
    assert results.loc[b"AAAAAA", "std_group1"] == 0.0  # No std with one sample
    assert results.loc[b"AAAAAA", "std_group2"] == 0.0  # No std with one sample


def test_single_sample_differential_abundance_identical():
    """Test single sample differential abundance calculation with identical values"""
    # Create identical DataFrames
    df = pd.DataFrame(
        {
            "sample1": [10, 20, 30],
        },
        index=[b"AAAAAA", b"TTTTTT", b"CCCCCC"],
    )

    results = SingleSampleDifferentialAbundance.calculate_differential_abundance_single(
        df, df
    )

    # All log2 fold changes should be 0
    assert all(abs(results["log2_fold_change"]) < 0.01)
    # All effect sizes should be 0
    assert all(abs(results["effect_size"]) < 0.01)


def test_single_sample_differential_abundance_empty():
    """Test single sample differential abundance calculation with empty DataFrame"""
    empty_df = pd.DataFrame()
    df = pd.DataFrame({"sample1": [1]}, index=[b"AAAAAA"])

    with pytest.raises(ValueError):
        SingleSampleDifferentialAbundance.calculate_differential_abundance_single(
            empty_df, df
        )
    with pytest.raises(ValueError):
        SingleSampleDifferentialAbundance.calculate_differential_abundance_single(
            df, empty_df
        )


def test_single_sample_differential_abundance_multiple_samples():
    """Test single sample differential abundance calculation with multiple samples"""
    df1 = pd.DataFrame(
        {
            "sample1": [10],
            "sample2": [20],
        },
        index=[b"AAAAAA"],
    )

    df2 = pd.DataFrame(
        {
            "sample3": [30],
        },
        index=[b"AAAAAA"],
    )

    with pytest.raises(ValueError):
        SingleSampleDifferentialAbundance.calculate_differential_abundance_single(
            df1, df2
        )


def test_single_sample_differential_abundance_with_allowed_list(partial_umi):
    """Test single sample differential abundance calculation with allowed list"""
    allowed_list = ["AAAAAA", "TTTTTT", "CCCCCC", "GGGGGG", "ATATAT"]
    result = SingleSampleDifferentialAbundance(
        [partial_umi, partial_umi], ["group1", "group2"], allowed_list=allowed_list
    ).run()
    assert result is not None
    assert isinstance(result, pd.DataFrame)
    assert all(col in result.columns for col in ["log2_fold_change", "effect_size"])
    assert "p_value" not in result.columns  # No p-values in single sample comparison
