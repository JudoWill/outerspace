"""Tests for UMI efficiency rate calculations"""

import pytest
from outerspace.stats import UMIEfficiencyRate


def test_efficiency_rate():
    """Test efficiency rate calculation with mixed allowed and banned UMIs"""
    counts = {
        b"AAAAAA": 10,  # Allowed
        b"TTTTTT": 20,  # Allowed
        b"CCCCCC": 30,  # Banned
        b"GGGGGG": 40,  # Banned
    }
    allowed_list = [b"AAAAAA", b"TTTTTT"]
    efficiency_rate = UMIEfficiencyRate.calculate_efficiency_rate(counts, allowed_list)
    # 30 allowed reads out of 100 total reads
    assert abs(efficiency_rate - 0.3) < 0.01


def test_efficiency_rate_all_allowed():
    """Test efficiency rate calculation with all UMIs allowed"""
    counts = {
        b"AAAAAA": 10,
        b"TTTTTT": 20,
    }
    allowed_list = [b"AAAAAA", b"TTTTTT"]
    efficiency_rate = UMIEfficiencyRate.calculate_efficiency_rate(counts, allowed_list)
    assert abs(efficiency_rate - 1.0) < 0.01  # Should be 1


def test_efficiency_rate_none_allowed():
    """Test efficiency rate calculation with no allowed UMIs"""
    counts = {
        b"CCCCCC": 30,
        b"GGGGGG": 40,
    }
    allowed_list = [b"AAAAAA", b"TTTTTT"]
    efficiency_rate = UMIEfficiencyRate.calculate_efficiency_rate(counts, allowed_list)
    assert abs(efficiency_rate) < 0.01  # Should be 0


def test_efficiency_rate_empty_counts():
    """Test efficiency rate calculation with empty counts"""
    counts = {}
    allowed_list = [b"AAAAAA", b"TTTTTT"]
    efficiency_rate = UMIEfficiencyRate.calculate_efficiency_rate(counts, allowed_list)
    assert efficiency_rate is None


def test_efficiency_rate_empty_allowed():
    """Test efficiency rate calculation with empty allowed list"""
    counts = {
        b"AAAAAA": 10,
        b"TTTTTT": 20,
    }
    allowed_list = []
    efficiency_rate = UMIEfficiencyRate.calculate_efficiency_rate(counts, allowed_list)
    assert abs(efficiency_rate) < 0.01  # Should be 0


def test_efficiency_rate_zero_counts():
    """Test efficiency rate calculation with zero counts"""
    counts = {
        b"AAAAAA": 0,
        b"TTTTTT": 0,
    }
    allowed_list = [b"AAAAAA", b"TTTTTT"]
    efficiency_rate = UMIEfficiencyRate.calculate_efficiency_rate(counts, allowed_list)
    assert efficiency_rate is None


def test_efficiency_rate_end_to_end(partial_umi):
    """Test end-to-end efficiency rate calculation"""
    allowed_list = ["AAAAAA", "TTTTTT", "CCCCCC", "GGGGGG", "ATATAT"]
    result = UMIEfficiencyRate.calculate(partial_umi, allowed_list=allowed_list)
    assert result is not None
    assert 0 <= result <= 1


def test_efficiency_rate_no_allowed_list(partial_umi):
    """Test efficiency rate calculation without allowed list"""
    result = UMIEfficiencyRate.calculate(partial_umi)
    assert result is None
