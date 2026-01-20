"""Tests for UMI alignment functionality using spoa.

This module contains unit tests for the alignment logic in outerspace.align,
testing the align_umi_sequences function independently of the CLI.
"""

import pytest
from outerspace.align import align_umi_sequences

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


def test_basic_alignment():
    """Test basic alignment of similar sequences"""
    sequences = ["AACTTATA", "AACTTATG", "AACTATA"]
    aligned = align_umi_sequences(sequences)

    # Should return same number of sequences
    assert len(aligned) == 3

    # All aligned sequences should have the same length
    aligned_lengths = set(len(seq) for seq in aligned)
    assert len(aligned_lengths) == 1, "All aligned sequences should have the same length"

    # Verify sequences are not empty
    assert all(len(seq) > 0 for seq in aligned)


def test_alignment_single_sequence():
    """Test alignment with a single sequence"""
    sequences = ["AACTTATA"]
    aligned = align_umi_sequences(sequences)

    # Single sequence should be returned as-is
    assert len(aligned) == 1
    assert aligned[0] == "AACTTATA"


def test_alignment_empty_input():
    """Test alignment with empty input"""
    sequences = []
    aligned = align_umi_sequences(sequences)

    # Should return empty list
    assert aligned == []


def test_alignment_empty_sequences():
    """Test alignment with sequences containing empty strings"""
    sequences = ["AACTTATA", "", "AACTATA", ""]
    aligned = align_umi_sequences(sequences)

    # Empty sequences should be filtered out
    assert len(aligned) == 2
    assert all(len(seq) > 0 for seq in aligned)
    assert all(len(seq) == len(aligned[0]) for seq in aligned)


def test_alignment_different_lengths():
    """Test alignment of sequences with different lengths"""
    sequences = ["AACTTATA", "AACTTATAG", "AACT"]
    aligned = align_umi_sequences(sequences)

    # Should align sequences of different lengths
    assert len(aligned) == 3
    # All should have the same length after alignment
    aligned_lengths = set(len(seq) for seq in aligned)
    assert len(aligned_lengths) == 1


def test_alignment_with_parameters():
    """Test alignment with custom parameters"""
    sequences = ["AACTTATA", "AACTTATG", "AACTATA"]

    # Test with different match score
    aligned1 = align_umi_sequences(sequences, match=10)
    assert len(aligned1) == 3
    assert all(len(seq) == len(aligned1[0]) for seq in aligned1)

    # Test with different mismatch penalty
    aligned2 = align_umi_sequences(sequences, mismatch=-2)
    assert len(aligned2) == 3
    assert all(len(seq) == len(aligned2[0]) for seq in aligned2)

    # Test with different gap penalty
    aligned3 = align_umi_sequences(sequences, gap=-5)
    assert len(aligned3) == 3
    assert all(len(seq) == len(aligned3[0]) for seq in aligned3)


def test_alignment_algorithm_local():
    """Test alignment with local algorithm (algorithm=0)"""
    sequences = ["AACTTATA", "AACTTATG", "AACTATA"]
    aligned = align_umi_sequences(sequences, algorithm=0)

    assert len(aligned) == 3
    assert all(len(seq) == len(aligned[0]) for seq in aligned)


def test_alignment_algorithm_global():
    """Test alignment with global algorithm (algorithm=1)"""
    sequences = ["AACTTATA", "AACTTATG", "AACTATA"]
    aligned = align_umi_sequences(sequences, algorithm=1)

    assert len(aligned) == 3
    assert all(len(seq) == len(aligned[0]) for seq in aligned)


def test_alignment_algorithm_semi_global():
    """Test alignment with semi-global algorithm (algorithm=2)"""
    sequences = ["AACTTATA", "AACTTATG", "AACTATA"]
    aligned = align_umi_sequences(sequences, algorithm=2)

    assert len(aligned) == 3
    assert all(len(seq) == len(aligned[0]) for seq in aligned)


def test_alignment_invalid_algorithm():
    """Test that invalid algorithm raises ValueError"""
    sequences = ["AACTTATA", "AACTTATG"]

    with pytest.raises(ValueError, match="Algorithm must be 0, 1, or 2"):
        align_umi_sequences(sequences, algorithm=3)

    with pytest.raises(ValueError, match="Algorithm must be 0, 1, or 2"):
        align_umi_sequences(sequences, algorithm=-1)


def test_alignment_identical_sequences():
    """Test alignment of identical sequences"""
    sequences = ["AACTTATA", "AACTTATA", "AACTTATA"]
    aligned = align_umi_sequences(sequences)

    assert len(aligned) == 3
    # All should be identical after alignment
    assert all(seq == aligned[0] for seq in aligned)


def test_alignment_very_different_sequences():
    """Test alignment of very different sequences"""
    sequences = ["AAAAA", "TTTTT", "CCCCC"]
    aligned = align_umi_sequences(sequences)

    assert len(aligned) == 3
    # Should still produce aligned sequences of same length
    assert all(len(seq) == len(aligned[0]) for seq in aligned)


def test_alignment_contains_gaps():
    """Test that alignment can introduce gaps"""
    sequences = ["AACTTATA", "AACTATA"]  # One sequence is shorter
    aligned = align_umi_sequences(sequences)

    assert len(aligned) == 2
    assert all(len(seq) == len(aligned[0]) for seq in aligned)
    # The shorter sequence should have gaps inserted
    assert len(aligned[0]) == len(aligned[1])


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
