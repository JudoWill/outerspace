"""Tests for Spearman correlation calculations"""

import pytest
import numpy as np
from outerspace.stats import SpearmanCorrelation


def test_spearman_calculation():
    """Test Spearman correlation calculation with known counts"""
    counts1 = {
        b"AAAAAA": 10,
        b"TTTTTT": 20,
        b"CCCCCC": 30,
    }
    counts2 = {
        b"AAAAAA": 20,
        b"TTTTTT": 40,
        b"CCCCCC": 60,
    }
    correlation = SpearmanCorrelation.calculate_spearman(counts1, counts2)
    # Perfect positive correlation (doubled values)
    assert abs(correlation - 1.0) < 0.01


def test_spearman_negative():
    """Test Spearman correlation calculation with negative correlation"""
    counts1 = {
        b"AAAAAA": 10,
        b"TTTTTT": 20,
        b"CCCCCC": 30,
    }
    counts2 = {
        b"AAAAAA": 30,
        b"TTTTTT": 20,
        b"CCCCCC": 10,
    }
    correlation = SpearmanCorrelation.calculate_spearman(counts1, counts2)
    # Perfect negative correlation
    assert abs(correlation + 1.0) < 0.01


def test_spearman_no_overlap():
    """Test Spearman correlation calculation with no overlapping UMIs"""
    counts1 = {
        b"AAAAAA": 10,
        b"TTTTTT": 20,
    }
    counts2 = {
        b"CCCCCC": 30,
        b"GGGGGG": 40,
    }
    correlation = SpearmanCorrelation.calculate_spearman(counts1, counts2)
    assert abs(correlation) < 0.01


def test_spearman_empty():
    """Test Spearman correlation calculation with empty counts"""
    counts1 = {}
    counts2 = {}
    correlation = SpearmanCorrelation.calculate_spearman(counts1, counts2)
    assert abs(correlation) < 0.01


def test_spearman_with_allowed_list(partial_umi):
    """Test Spearman correlation calculation with allowed list"""
    allowed_list = ["AAAAAA", "TTTTTT", "CCCCCC", "GGGGGG", "ATATAT"]
    result = SpearmanCorrelation.calculate(
        partial_umi, partial_umi, allowed_list=allowed_list
    )
    assert result is not None
    assert -1 <= result <= 1
