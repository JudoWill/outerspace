"""Tests for UMI redundancy calculations"""

import pytest
from outerspace.stats import UMIRedundancy


def test_redundancy_calculation():
    """Test redundancy calculation with known counts"""
    counts = {
        b"AAAAAA": 10,
        b"TTTTTT": 20,
        b"CCCCCC": 30,
    }
    redundancy = UMIRedundancy.calculate_redundancy(counts)
    # Total reads = 60, unique UMIs = 3
    # Redundancy = 60/3 = 20
    assert abs(redundancy - 20.0) < 0.01


def test_redundancy_perfect():
    """Test redundancy calculation with perfect redundancy (1 read per UMI)"""
    counts = {
        b"AAAAAA": 1,
        b"TTTTTT": 1,
        b"CCCCCC": 1,
    }
    redundancy = UMIRedundancy.calculate_redundancy(counts)
    assert abs(redundancy - 1.0) < 0.01


def test_redundancy_empty_counts():
    """Test redundancy calculation with empty counts"""
    counts = {}
    redundancy = UMIRedundancy.calculate_redundancy(counts)
    assert redundancy is None


def test_redundancy_zero_counts():
    """Test redundancy calculation with zero counts"""
    counts = {
        b"AAAAAA": 0,
        b"TTTTTT": 0,
        b"CCCCCC": 0,
    }
    redundancy = UMIRedundancy.calculate_redundancy(counts)
    assert 0 <= redundancy < 0.01


def test_redundancy_mixed_zero_counts():
    """Test redundancy calculation with some zero counts"""
    counts = {
        b"AAAAAA": 10,
        b"TTTTTT": 0,
        b"CCCCCC": 20,
    }
    redundancy = UMIRedundancy.calculate_redundancy(counts)
    # Total reads = 30, unique UMIs = 3
    # Redundancy = 30/3 = 10
    assert abs(redundancy - 10.0) < 0.01


def test_redundancy_end_to_end(equal_umi):
    """Test end-to-end redundancy calculation"""
    result = UMIRedundancy.calculate(equal_umi)
    assert result is not None
    assert result >= 1.0


def test_redundancy_with_allowed_list(partial_umi):
    """Test redundancy calculation with allowed list"""
    allowed_list = ["AAAAAA", "TTTTTT", "CCCCCC", "GGGGGG", "ATATAT"]
    result = UMIRedundancy.calculate(partial_umi, allowed_list=allowed_list)
    assert result is not None
    assert result >= 1.0

