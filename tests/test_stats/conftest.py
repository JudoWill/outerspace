"""Shared fixtures for UMI statistics tests"""

import pytest
from outerspace.umi import UMI


@pytest.fixture
def empty_umi():
    """Empty UMI object with mapping created"""
    umi = UMI(mismatches=0)
    umi.create_mapping()
    return umi


@pytest.fixture
def equal_umi():
    """UMI object with equal counts for all sequences"""
    umi = UMI(mismatches=0)
    sequences = [
        "AAAAAA",
        "TTTTTT",
        "CCCCCC",
        "GGGGGG",
        "ATATAT",
        "TATATA",
        "CGCGCG",
        "GCGCGC",
        "ACACAC",
        "TGTGTG",
    ]
    for seq in sequences:
        umi.consume(seq)
    umi.create_mapping()
    return umi


@pytest.fixture
def unequal_umi():
    """UMI object with highly unequal distribution"""
    umi = UMI(mismatches=0)
    sequences = [
        "AAAAAA",
        "TTTTTT",
        "CCCCCC",
        "GGGGGG",
        "ATATAT",
        "TATATA",
        "CGCGCG",
        "GCGCGC",
        "ACACAC",
        "TGTGTG",
    ]
    # Add 100 counts of the first sequence
    for _ in range(100):
        umi.consume(sequences[0])
    # Add 1 count of each other sequence
    for seq in sequences[1:]:
        umi.consume(seq)
    umi.create_mapping()
    return umi


@pytest.fixture
def moderate_umi():
    """UMI object with moderate inequality"""
    umi = UMI(mismatches=0)
    sequences = [
        "AAAAAA",
        "TTTTTT",
        "CCCCCC",
        "GGGGGG",
        "ATATAT",
        "TATATA",
        "CGCGCG",
        "GCGCGC",
        "ACACAC",
        "TGTGTG",
    ]
    counts = [1, 2, 3, 4, 5, 5, 4, 3, 2, 1]
    for seq, count in zip(sequences, counts):
        for _ in range(count):
            umi.consume(seq)
    umi.create_mapping()
    return umi


@pytest.fixture
def merged_umi():
    """UMI object with sequences that will be merged"""
    umi = UMI(mismatches=1)
    umi.consume("ATCG", 10)
    umi.consume("ATCC", 5)  # Will be merged with ATCG
    umi.consume("GCTA", 3)
    umi.create_mapping()
    return umi


@pytest.fixture
def partial_umi():
    """UMI object with partial distribution and some missing keys"""
    umi = UMI(mismatches=0)
    sequences = ["AAAAAA", "TTTTTT", "CCCCCC", "GGGGGG", "ATATAT"]
    counts = [10, 8, 6, 4, 2]  # Unequal distribution
    for seq, count in zip(sequences, counts):
        for _ in range(count):
            umi.consume(seq)
    umi.create_mapping()
    return umi
