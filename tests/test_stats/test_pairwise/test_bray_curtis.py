"""Tests for Bray-Curtis dissimilarity calculations"""

import pytest
import numpy as np
from outerspace.stats import BrayCurtisDissimilarity

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

def test_bray_curtis_calculation():
    """Test Bray-Curtis dissimilarity calculation with known counts"""
    counts1 = {
        b"AAAAAA": 10,
        b"TTTTTT": 20,
        b"CCCCCC": 30,
    }
    counts2 = {
        b"AAAAAA": 20,
        b"TTTTTT": 10,
        b"GGGGGG": 30,
    }
    dissimilarity = BrayCurtisDissimilarity.calculate_bray_curtis(counts1, counts2)
    # Total difference = |10-20| + |20-10| + |30-0| + |0-30| = 80
    # Total sum = (10+20) + (20+10) + (30+0) + (0+30) = 120
    # Dissimilarity = 80/120 = 0.666...

    expected_dissimilarity = 0.6666
    assert abs(dissimilarity - expected_dissimilarity) < 0.01


def test_bray_curtis_identical():
    """Test Bray-Curtis dissimilarity calculation with identical counts"""
    counts1 = {
        b"AAAAAA": 10,
        b"TTTTTT": 20,
        b"CCCCCC": 30,
    }
    counts2 = {
        b"AAAAAA": 10,
        b"TTTTTT": 20,
        b"CCCCCC": 30,
    }
    dissimilarity = BrayCurtisDissimilarity.calculate_bray_curtis(counts1, counts2)
    assert abs(dissimilarity) < 0.01


def test_bray_curtis_empty():
    """Test Bray-Curtis dissimilarity calculation with empty counts"""
    counts1 = {}
    counts2 = {}
    dissimilarity = BrayCurtisDissimilarity.calculate_bray_curtis(counts1, counts2)
    assert abs(dissimilarity) < 0.01


def test_bray_curtis_zero_counts():
    """Test Bray-Curtis dissimilarity calculation with zero counts"""
    counts1 = {
        b"AAAAAA": 0,
        b"TTTTTT": 0,
    }
    counts2 = {
        b"AAAAAA": 0,
        b"TTTTTT": 0,
    }
    dissimilarity = BrayCurtisDissimilarity.calculate_bray_curtis(counts1, counts2)
    assert abs(dissimilarity) < 0.01


def test_bray_curtis_with_allowed_list(partial_umi):
    """Test Bray-Curtis dissimilarity calculation with allowed list"""
    allowed_list = ["AAAAAA", "TTTTTT", "CCCCCC", "GGGGGG", "ATATAT"]
    result = BrayCurtisDissimilarity.calculate(
        partial_umi, partial_umi, allowed_list=allowed_list
    )
    assert result is not None
    assert 0 <= result <= 1


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.