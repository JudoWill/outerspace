"""Tests for paired t-test differential abundance calculations"""

import pytest
import numpy as np
import pandas as pd
from outerspace.stats.differential import PairedTTestDifferentialAbundance


def test_paired_ttest_differential_abundance():
    """Test paired t-test differential abundance calculation with known counts"""
    # Create test DataFrames for pre/post samples
    pre_df = pd.DataFrame(
        {
            "sample1": [10, 20, 30],
            "sample2": [15, 25, 35],
        },
        index=[b"AAAAAA", b"TTTTTT", b"CCCCCC"],
    )

    post_df = pd.DataFrame(
        {
            "sample1": [20, 30, 40],  # All values increased by 10
            "sample2": [25, 35, 45],
        },
        index=[b"AAAAAA", b"TTTTTT", b"CCCCCC"],
    )

    results = PairedTTestDifferentialAbundance.calculate_differential_abundance_paired(
        pre_df, post_df
    )

    # Check results
    assert all(results["log2_fold_change"] > 0)  # All should be positive
    assert all(results["effect_size"] > 0)  # All should be positive
    assert all((0 <= results["p_value"]) & (results["p_value"] <= 1))  # Valid p-values
    assert all(
        results["mean_post"] > results["mean_pre"]
    )  # Post means should be higher


def test_paired_ttest_differential_abundance_varying_differences():
    """Test paired t-test differential abundance calculation with varying differences between genes"""
    # Create test DataFrames for pre/post samples with different fold changes
    pre_df = pd.DataFrame(
        {
            "sample1": [10, 20, 30, 40],
            "sample2": [15, 25, 35, 45],
        },
        index=[b"AAAAAA", b"TTTTTT", b"CCCCCC", b"GGGGGG"],
    )

    post_df = pd.DataFrame(
        {
            "sample1": [30, 40, 45, 60],  # Different increases: 20, 20, 15, 20
            "sample2": [35, 45, 50, 65],
        },
        index=[b"AAAAAA", b"TTTTTT", b"CCCCCC", b"GGGGGG"],
    )

    results = PairedTTestDifferentialAbundance.calculate_differential_abundance_paired(
        pre_df, post_df
    )

    # Check results
    assert all(results["log2_fold_change"] > 0)  # All should be positive
    assert all(results["effect_size"] > 0)  # All should be positive
    assert all((0 <= results["p_value"]) & (results["p_value"] <= 1))  # Valid p-values
    assert all(
        results["mean_post"] > results["mean_pre"]
    )  # Post means should be higher

    # Check that effect sizes are different between genes
    effect_sizes = results["effect_size"].values
    assert len(set(effect_sizes)) > 1  # Should have different effect sizes

    # Check that the gene with smallest increase (CCCCCC) has smallest effect size
    min_effect_idx = np.argmin(effect_sizes)
    assert results.index[min_effect_idx] == b"CCCCCC"

    # Check that the gene with largest increase (AAAAAA) has largest effect size
    max_effect_idx = np.argmax(effect_sizes)
    assert results.index[max_effect_idx] == b"AAAAAA"


def test_paired_ttest_differential_abundance_identical():
    """Test paired t-test differential abundance calculation with identical values"""
    # Create identical DataFrames
    df = pd.DataFrame(
        {
            "sample1": [10, 20, 30],
            "sample2": [10, 20, 30],
        },
        index=[b"AAAAAA", b"TTTTTT", b"CCCCCC"],
    )

    results = PairedTTestDifferentialAbundance.calculate_differential_abundance_paired(
        df, df
    )

    # All log2 fold changes should be 0
    assert all(abs(results["log2_fold_change"]) < 0.01)
    # All effect sizes should be 0
    assert all(abs(results["effect_size"]) < 0.01)
    # All p-values should be 1 (identical distributions)
    assert all(abs(results["p_value"] - 1.0) < 0.01)


def test_paired_ttest_differential_abundance_different():
    """Test paired t-test differential abundance calculation with different values"""
    # Create DataFrames with different values
    pre_df = pd.DataFrame(
        {
            "sample1": [10, 20, 30],
            "sample2": [15, 25, 35],
        },
        index=[b"AAAAAA", b"TTTTTT", b"CCCCCC"],
    )

    post_df = pd.DataFrame(
        {
            "sample1": [
                30,
                20,
                10,
            ],  # Values are different but not consistently higher/lower
            "sample2": [35, 25, 15],
        },
        index=[b"AAAAAA", b"TTTTTT", b"CCCCCC"],
    )

    results = PairedTTestDifferentialAbundance.calculate_differential_abundance_paired(
        pre_df, post_df
    )

    # Check that p-values are calculated and valid
    assert all((0 <= results["p_value"]) & (results["p_value"] <= 1))  # Valid p-values
    assert not all(abs(results["p_value"] - 1.0) < 0.01)  # Not all p-values should be 1

    # Check that effect sizes reflect the differences
    assert not all(results["effect_size"] == 0)  # Not all effect sizes should be 0
    assert not all(
        results["effect_size"] > 0
    )  # Not all effect sizes should be positive
    assert not all(
        results["effect_size"] < 0
    )  # Not all effect sizes should be negative


def test_paired_ttest_differential_abundance_empty():
    """Test paired t-test differential abundance calculation with empty DataFrame"""
    empty_df = pd.DataFrame()
    df = pd.DataFrame({"sample1": [1, 2, 3]}, index=[b"AAAAAA", b"TTTTTT", b"CCCCCC"])

    with pytest.raises(ValueError):
        PairedTTestDifferentialAbundance.calculate_differential_abundance_paired(
            empty_df, df
        )
    with pytest.raises(ValueError):
        PairedTTestDifferentialAbundance.calculate_differential_abundance_paired(
            df, empty_df
        )


def test_paired_ttest_differential_abundance_mismatched():
    """Test paired t-test differential abundance calculation with mismatched samples"""
    pre_df = pd.DataFrame(
        {
            "sample1": [10, 20, 30],
            "sample2": [15, 25, 35],
        },
        index=[b"AAAAAA", b"TTTTTT", b"CCCCCC"],
    )

    post_df = pd.DataFrame(
        {
            "sample1": [20, 30, 40],
        },
        index=[b"AAAAAA", b"TTTTTT", b"CCCCCC"],
    )

    with pytest.raises(ValueError):
        PairedTTestDifferentialAbundance.calculate_differential_abundance_paired(
            pre_df, post_df
        )


def test_paired_ttest_differential_abundance_with_allowed_list(partial_umi):
    """Test paired t-test differential abundance calculation with allowed list"""
    allowed_list = ["AAAAAA", "TTTTTT", "CCCCCC", "GGGGGG", "ATATAT"]
    result = PairedTTestDifferentialAbundance(
        [partial_umi, partial_umi], ["group1", "group2"], allowed_list=allowed_list
    ).run()
    assert result is not None
    assert isinstance(result, pd.DataFrame)
    assert all(
        col in result.columns for col in ["log2_fold_change", "effect_size", "p_value"]
    )
