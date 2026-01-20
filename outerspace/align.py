"""UMI sequence alignment functionality using spoa.

This module provides functions for aligning UMI sequences using the spoa
library. It supports configurable alignment parameters including match scores,
mismatch penalties, gap penalties, and different alignment algorithms.
"""

import logging
from typing import List

try:
    from spoa import poa
except ImportError:
    raise ImportError(
        "spoa is required for alignment functionality. "
        "Install it with: pip install spoa"
    )

# Set up logging
logger = logging.getLogger(__name__)

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


def align_umi_sequences(
    sequences: List[str],
    match: int = 5,
    mismatch: int = -4,
    gap: int = -8,
    algorithm: int = 1,
) -> List[str]:
    """Align UMI sequences using spoa (Partial Order Alignment).

    This function performs multiple sequence alignment on a list of UMI sequences
    using the spoa library. It returns aligned sequences where all sequences have
    the same length (with gaps inserted as needed).

    Parameters
    ----------
    sequences : List[str]
        List of UMI sequences to align
    match : int, default=5
        Score for matching bases
    mismatch : int, default=-4
        Penalty for mismatching bases (should be negative)
    gap : int, default=-8
        Penalty for gaps/indels (should be negative)
    algorithm : int, default=1
        Alignment algorithm type:
        - 0: Local alignment (Smith-Waterman)
        - 1: Global alignment (Needleman-Wunsch)
        - 2: Semi-global alignment (overlap)

    Returns
    -------
    List[str]
        List of aligned sequences (all with the same length)

    Raises
    ------
    ValueError
        If sequences list is empty or contains invalid sequences
    ImportError
        If spoa is not installed

    Examples
    --------
    >>> sequences = ["AACTTATA", "AACTTATG", "AACTATA"]
    >>> aligned = align_umi_sequences(sequences)
    >>> len(aligned) == 3
    True
    >>> all(len(seq) == len(aligned[0]) for seq in aligned)
    True
    """
    if not sequences:
        logger.warning("Empty sequence list provided")
        return []

    # Filter out empty sequences
    non_empty_sequences = [seq for seq in sequences if seq]
    if not non_empty_sequences:
        logger.warning("All sequences are empty")
        return []

    # If only one sequence, return it as-is (no alignment needed)
    if len(non_empty_sequences) == 1:
        logger.debug("Only one sequence provided, returning as-is")
        return non_empty_sequences

    # Validate algorithm parameter
    if algorithm not in [0, 1, 2]:
        raise ValueError(f"Algorithm must be 0, 1, or 2, got {algorithm}")

    logger.debug(
        f"Aligning {len(non_empty_sequences)} sequences with "
        f"match={match}, mismatch={mismatch}, gap={gap}, algorithm={algorithm}"
    )

    try:
        # spoa.poa returns (consensus, msa) tuple
        # We use the msa (multiple sequence alignment) which is a list of aligned sequences
        # Parameter names: m=match, n=mismatch, g=gap
        # genmsa=True ensures we get the MSA output
        consensus, msa = poa(
            non_empty_sequences,
            algorithm=algorithm,
            genmsa=True,
            m=match,
            n=mismatch,
            g=gap,
        )

        if not msa:
            logger.warning("spoa returned empty alignment")
            return non_empty_sequences

        # Verify all aligned sequences have the same length
        if msa:
            aligned_length = len(msa[0])
            if not all(len(seq) == aligned_length for seq in msa):
                logger.warning(
                    "Aligned sequences have different lengths, "
                    "this should not happen with spoa"
                )

        logger.debug(f"Successfully aligned {len(msa)} sequences to length {len(msa[0])}")
        return msa

    except Exception as e:
        logger.error(f"Error during alignment: {e}")
        raise ValueError(f"Alignment failed: {e}") from e


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
