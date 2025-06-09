"""Tests for UMI recovery rate calculations"""

import pytest
import numpy as np
from outerspace.stats import UMIRecoveryRate


def test_recovery_rate_limited():
    """Test recovery rate calculation with allowed list"""
    counts = {
        b"AAAAAA": 10,
        b"TTTTTT": 20,
        b"CCCCCC": 30,
    }
    allowed_list = [b"AAAAAA", b"TTTTTT", b"GGGGGG", b"ATATAT"]
    recovery_rate = UMIRecoveryRate.calculate_recovery_rate_limited(counts, allowed_list)
    # 2 out of 4 allowed UMIs are present
    assert abs(recovery_rate - 0.5) < 0.01


def test_recovery_rate_limited_empty():
    """Test recovery rate calculation with empty counts"""
    counts = {}
    allowed_list = [b"AAAAAA", b"TTTTTT"]
    recovery_rate = UMIRecoveryRate.calculate_recovery_rate_limited(counts, allowed_list)
    assert abs(recovery_rate) < 0.01  # Should be 0


def test_recovery_rate_limited_all_present():
    """Test recovery rate calculation with all allowed UMIs present"""
    counts = {
        b"AAAAAA": 10,
        b"TTTTTT": 20,
    }
    allowed_list = [b"AAAAAA", b"TTTTTT"]
    recovery_rate = UMIRecoveryRate.calculate_recovery_rate_limited(counts, allowed_list)
    assert abs(recovery_rate - 1.0) < 0.01  # Should be 1



def test_recovery_rate_with_allowed_list(partial_umi):
    """Test recovery rate calculation with allowed list"""
    allowed_list = ["AAAAAA", "TTTTTT", "CCCCCC", "GGGGGG", "ATATAT"]
    result = UMIRecoveryRate.calculate(partial_umi, allowed_list=allowed_list)
    assert 'umi_recovery_rate' in result
    assert result['umi_recovery_rate'] is not None
    assert result['umi_recovery_rate'] == 1 