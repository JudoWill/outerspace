"""Tests for UMI error rate calculations"""

import pytest
from outerspace.stats import UMIErrorRate


def test_hamming_distance():
    """Test Hamming distance calculation between sequences"""
    seq1 = b"AAAAAA"
    seq2 = b"ATAAAA"
    distance = UMIErrorRate.hamming_distance(seq1, seq2)
    assert distance == 1

    seq1 = b"AAAAAA"
    seq2 = b"TTTTTT"
    distance = UMIErrorRate.hamming_distance(seq1, seq2)
    assert distance == 6

    seq1 = b"AAAAAA"
    seq2 = b"AAAAAA"
    distance = UMIErrorRate.hamming_distance(seq1, seq2)
    assert distance == 0


def test_hamming_distance_different_lengths():
    """Test Hamming distance calculation with sequences of different lengths"""
    seq1 = b"AAAAAA"
    seq2 = b"AAAAA"
    with pytest.raises(ValueError):
        UMIErrorRate.hamming_distance(seq1, seq2)


def test_error_rate_calculation():
    """Test error rate calculation with known mismatches"""
    mapping = {
        b"AAAAAA": b"AAAAAA",  # No error
        b"ATAAAA": b"AAAAAA",  # 1 mismatch
        b"TTAAAA": b"AAAAAA",  # 2 mismatches
    }
    counts = {
        b"AAAAAA": 10,
        b"ATAAAA": 20,
        b"TTAAAA": 30,
    }
    error_rate = UMIErrorRate.calculate_error_rate(mapping, counts)
    # Total mismatches = (1 * 20) + (2 * 30) = 80
    # Total reads = 60
    # Error rate = 80/60 = 1.33
    assert abs(error_rate - 1.33) < 0.01


def test_error_rate_no_errors():
    """Test error rate calculation with no errors"""
    mapping = {
        b"AAAAAA": b"AAAAAA",
        b"TTTTTT": b"TTTTTT",
    }
    counts = {
        b"AAAAAA": 10,
        b"TTTTTT": 20,
    }
    error_rate = UMIErrorRate.calculate_error_rate(mapping, counts)
    assert abs(error_rate) < 0.01  # Should be 0


def test_error_rate_empty_mapping():
    """Test error rate calculation with empty mapping"""
    mapping = {}
    counts = {
        b"AAAAAA": 10,
        b"TTTTTT": 20,
    }
    error_rate = UMIErrorRate.calculate_error_rate(mapping, counts)
    assert error_rate is None


def test_error_rate_empty_counts():
    """Test error rate calculation with empty counts"""
    mapping = {
        b"AAAAAA": b"AAAAAA",
        b"TTTTTT": b"TTTTTT",
    }
    counts = {}
    error_rate = UMIErrorRate.calculate_error_rate(mapping, counts)
    assert error_rate is None


def test_error_rate_end_to_end(merged_umi):
    """Test end-to-end error rate calculation"""
    result = UMIErrorRate.calculate(merged_umi)
    assert result is not None
    assert result >= 0


def test_error_rate_no_mapping(empty_umi):
    """Test error rate calculation with no mapping"""
    result = UMIErrorRate.calculate(empty_umi)
    assert result is None 