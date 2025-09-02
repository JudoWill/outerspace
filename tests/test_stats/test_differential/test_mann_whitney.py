"""Tests for Mann-Whitney differential abundance calculations"""

import pytest
import numpy as np
import pandas as pd
from outerspace.stats.differential import MannWhitneyDifferentialAbundance

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


def test_mann_whitney_differential_abundance():
    """Test Mann-Whitney differential abundance calculation with known counts"""
    # Create test DataFrames for each group
    group1_df = pd.DataFrame(
        {
            "sample1": [10, 20, 30, 0],
            "sample2": [15, 25, 35, 5],
        },
        index=[b"AAAAAA", b"TTTTTT", b"CCCCCC", b"GGGGGG"],
    )

    group2_df = pd.DataFrame(
        {
            "sample3": [20, 10, 0, 30],
            "sample4": [25, 15, 5, 35],
        },
        index=[b"AAAAAA", b"TTTTTT", b"CCCCCC", b"GGGGGG"],
    )

    results = MannWhitneyDifferentialAbundance.calculate_differential_abundance(
        group1_df, group2_df
    )

    mean_group1 = (10 + 15) / 2
    mean_group2 = (20 + 25) / 2

    fc = mean_group2 / mean_group1
    log2_fc = np.log2(fc)

    # Check results for AAAAAA
    assert abs(results.loc[b"AAAAAA", "log2_fold_change"] - log2_fc) < 0.01
    assert abs(results.loc[b"AAAAAA", "mean_group1"] - mean_group1) < 0.01
    assert abs(results.loc[b"AAAAAA", "mean_group2"] - mean_group2) < 0.01
    assert all((0 <= results["p_value"]) & (results["p_value"] <= 1))  # Valid p-values
    assert results.loc[b"AAAAAA", "effect_size"] > 0


def test_mann_whitney_differential_abundance_identical():
    """Test Mann-Whitney differential abundance calculation with identical groups"""
    # Create identical DataFrames
    df = pd.DataFrame(
        {
            "sample1": [10, 20, 30],
            "sample2": [10, 20, 30],
        },
        index=[b"AAAAAA", b"TTTTTT", b"CCCCCC"],
    )

    results = MannWhitneyDifferentialAbundance.calculate_differential_abundance(df, df)

    # All log2 fold changes should be 0
    assert all(abs(results["log2_fold_change"]) < 0.01)
    # All effect sizes should be 0
    assert all(abs(results["effect_size"]) < 0.01)
    # All p-values should be 1 (identical distributions)
    assert all(abs(results["p_value"] - 1.0) < 0.01)


def test_mann_whitney_differential_abundance_empty():
    """Test Mann-Whitney differential abundance calculation with empty DataFrame"""
    empty_df = pd.DataFrame()
    df = pd.DataFrame({"sample1": [1, 2, 3]}, index=[b"AAAAAA", b"TTTTTT", b"CCCCCC"])

    with pytest.raises(ValueError):
        MannWhitneyDifferentialAbundance.calculate_differential_abundance(empty_df, df)
    with pytest.raises(ValueError):
        MannWhitneyDifferentialAbundance.calculate_differential_abundance(df, empty_df)


def test_mann_whitney_differential_abundance_zero_counts():
    """Test Mann-Whitney differential abundance calculation with zero counts"""
    zero_df = pd.DataFrame(
        {
            "sample1": [0, 0, 0],
            "sample2": [0, 0, 0],
        },
        index=[b"AAAAAA", b"TTTTTT", b"CCCCCC"],
    )

    results = MannWhitneyDifferentialAbundance.calculate_differential_abundance(
        zero_df, zero_df
    )
    assert all(results["log2_fold_change"] == 0)
    assert all(results["effect_size"] == 0)
    assert all(abs(results["p_value"] - 1.0) < 0.01)  # Identical zero distributions


def test_mann_whitney_differential_abundance_with_allowed_list(partial_umi):
    """Test Mann-Whitney differential abundance calculation with allowed list"""
    allowed_list = ["AAAAAA", "TTTTTT", "CCCCCC", "GGGGGG", "ATATAT"]
    result = MannWhitneyDifferentialAbundance(
        [partial_umi, partial_umi], ["group1", "group2"], allowed_list=allowed_list
    ).run()
    assert result is not None
    assert isinstance(result, pd.DataFrame)
    assert all(
        col in result.columns for col in ["log2_fold_change", "effect_size", "p_value"]
    )

    
    # Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.