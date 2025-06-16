"""Tests for Jaccard similarity calculations"""

import pytest
import numpy as np
from outerspace.stats import JaccardSimilarity


def test_jaccard_calculation():
    """Test Jaccard similarity calculation with known sets"""
    set1 = {b"AAAAAA", b"TTTTTT", b"CCCCCC"}
    set2 = {b"AAAAAA", b"TTTTTT", b"GGGGGG"}
    similarity = JaccardSimilarity.calculate_jaccard(set1, set2)
    # 2 shared UMIs out of 4 total unique UMIs
    assert abs(similarity - 0.5) < 0.01


def test_jaccard_perfect():
    """Test Jaccard similarity calculation with identical sets"""
    set1 = {b"AAAAAA", b"TTTTTT", b"CCCCCC"}
    set2 = {b"AAAAAA", b"TTTTTT", b"CCCCCC"}
    similarity = JaccardSimilarity.calculate_jaccard(set1, set2)
    assert abs(similarity - 1.0) < 0.01


def test_jaccard_no_overlap():
    """Test Jaccard similarity calculation with no overlap"""
    set1 = {b"AAAAAA", b"TTTTTT"}
    set2 = {b"CCCCCC", b"GGGGGG"}
    similarity = JaccardSimilarity.calculate_jaccard(set1, set2)
    assert abs(similarity) < 0.01


def test_jaccard_empty_sets():
    """Test Jaccard similarity calculation with empty sets"""
    set1 = set()
    set2 = set()
    similarity = JaccardSimilarity.calculate_jaccard(set1, set2)
    assert abs(similarity) < 0.01


def test_jaccard_with_allowed_list(partial_umi):
    """Test Jaccard similarity calculation with allowed list"""
    allowed_list = ["AAAAAA", "TTTTTT", "CCCCCC", "GGGGGG", "ATATAT"]
    result = JaccardSimilarity.calculate(partial_umi, partial_umi, allowed_list=allowed_list)
    assert result is not None
    assert 0 <= result <= 1 