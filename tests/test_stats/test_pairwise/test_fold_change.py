"""Tests for fold change calculations"""

import pytest
import numpy as np
from outerspace.stats import FoldChange


def test_fold_change_calculation():
    """Test fold change calculation with known counts"""
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
    fold_changes = FoldChange.calculate_fold_change(counts1, counts2)
    assert abs(fold_changes[b"AAAAAA"] - 2.0) < 0.01  # 20/10
    assert abs(fold_changes[b"TTTTTT"] - 0.5) < 0.01  # 10/20
    assert abs(fold_changes[b"CCCCCC"] - 0.0) < 0.01  # 0/30
    assert fold_changes[b"GGGGGG"] == float("inf")  # 30/0


def test_fold_change_identical():
    """Test fold change calculation with identical counts"""
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
    fold_changes = FoldChange.calculate_fold_change(counts1, counts2)
    assert all(abs(fc - 1.0) < 0.01 for fc in fold_changes.values())


def test_fold_change_empty():
    """Test fold change calculation with empty counts"""
    counts1 = {}
    counts2 = {}
    fold_changes = FoldChange.calculate_fold_change(counts1, counts2)
    assert len(fold_changes) == 0


def test_fold_change_zero_counts():
    """Test fold change calculation with zero counts"""
    counts1 = {
        b"AAAAAA": 0,
        b"TTTTTT": 0,
    }
    counts2 = {
        b"AAAAAA": 0,
        b"TTTTTT": 0,
    }
    fold_changes = FoldChange.calculate_fold_change(counts1, counts2)
    assert all(fc == 0.0 for fc in fold_changes.values())


def test_fold_change_with_allowed_list(partial_umi):
    """Test fold change calculation with allowed list"""
    allowed_list = ["AAAAAA", "TTTTTT", "CCCCCC", "GGGGGG", "ATATAT"]
    result = FoldChange.calculate(partial_umi, partial_umi, allowed_list=allowed_list)
    assert result is not None
    assert all(fc >= 0 for fc in result.values())
