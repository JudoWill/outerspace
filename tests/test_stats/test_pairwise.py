"""Tests for pairwise UMI library statistics"""

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
from outerspace.umi import UMI, UmiCollection
from outerspace.stats.multi import (
    JaccardSimilarity,
    BrayCurtisDissimilarity,
    FoldChange,
    SpearmanCorrelation,
)

# List of all pairwise stat classes to test
PAIRWISE_STAT_CLASSES = [
    JaccardSimilarity,
    BrayCurtisDissimilarity,
    FoldChange,
    SpearmanCorrelation,
]


def test_pairwise_stat_interface(empty_umi, equal_umi):
    """Test that pairwise stats accept UMI objects and return correct types"""
    for stat_class in PAIRWISE_STAT_CLASSES:
        # Test calculate method
        result = stat_class.calculate(empty_umi, equal_umi)
        assert isinstance(
            result, (float, dict)
        ), f"{stat_class.__name__}.calculate() should return float or dict"


def test_pairwise_stat_collection_interface(empty_umi, equal_umi):
    """Test that pairwise stats accept UmiCollection objects and return correct shapes"""
    # Create a collection with two UMIs
    collection = UmiCollection({"sample1": empty_umi, "sample2": equal_umi})

    for stat_class in PAIRWISE_STAT_CLASSES:
        # Test calculate_collection method
        results = stat_class.calculate_collection(collection)
        assert isinstance(
            results, dict
        ), f"{stat_class.__name__}.calculate_collection() should return a dict"
        assert (
            "sample1",
            "sample2",
        ) in results, (
            f"{stat_class.__name__}.calculate_collection() should include all pairs"
        )
        assert isinstance(
            results[("sample1", "sample2")], (float, dict)
        ), f"{stat_class.__name__}.calculate_collection() values should be float or dict"


def test_pairwise_stat_collection_multiple_samples(empty_umi, equal_umi, unequal_umi):
    """Test that pairwise stats handle multiple samples correctly"""
    # Create a collection with multiple UMIs
    collection = UmiCollection(
        {"sample1": empty_umi, "sample2": equal_umi, "sample3": unequal_umi}
    )

    for stat_class in PAIRWISE_STAT_CLASSES:
        # Test calculate_collection method
        results = stat_class.calculate_collection(collection)
        assert isinstance(
            results, dict
        ), f"{stat_class.__name__}.calculate_collection() should return a dict"

        # Check that we have all possible pairs
        expected_pairs = {
            ("sample1", "sample2"),
            ("sample1", "sample3"),
            ("sample2", "sample3"),
        }
        assert (
            set(results.keys()) == expected_pairs
        ), f"{stat_class.__name__}.calculate_collection() should include all pairs"

        # Check value types
        for value in results.values():
            assert isinstance(
                value, (float, dict)
            ), f"{stat_class.__name__}.calculate_collection() values should be float or dict"


def test_pairwise_stat_with_allowed_list(empty_umi, equal_umi):
    """Test that pairwise stats handle allowed_list parameter correctly"""
    allowed_list = ["AAAAAA", "TTTTTT"]

    for stat_class in PAIRWISE_STAT_CLASSES:
        # Test calculate method with allowed_list
        result = stat_class.calculate(empty_umi, equal_umi, allowed_list=allowed_list)
        assert isinstance(
            result, (float, dict)
        ), f"{stat_class.__name__}.calculate() should handle allowed_list parameter"

        # Test calculate_collection method with allowed_list
        collection = UmiCollection({"sample1": empty_umi, "sample2": equal_umi})
        results = stat_class.calculate_collection(collection, allowed_list=allowed_list)
        assert isinstance(
            results, dict
        ), f"{stat_class.__name__}.calculate_collection() should handle allowed_list parameter"
        assert (
            "sample1",
            "sample2",
        ) in results, (
            f"{stat_class.__name__}.calculate_collection() should include all pairs"
        )

        
# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.