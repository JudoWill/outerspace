"""Utility functions for UMI statistics calculations.

This module provides utility functions for common operations in UMI statistics
calculations, such as filtering counts by allowed lists and data preprocessing.
"""

import logging
from typing import Dict, List, Tuple

# Set up logging
logger = logging.getLogger(__name__)

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


def split_counts_by_allowed_list(
    counts: Dict[bytes, int], allowed_list: List[bytes], add_missing: bool = True
) -> Tuple[Dict[bytes, int], Dict[bytes, int], List[bytes]]:
    """Split counts by allowed list.

    This function separates UMI counts into allowed and banned categories based
    on a provided allowed list, and optionally adds missing UMIs with zero counts.

    Parameters
    ----------
    counts : Dict[bytes, int]
        Dictionary of UMI counts
    allowed_list : List[bytes]
        List of allowed UMIs
    add_missing : bool, default=True
        Whether to add missing UMIs to the allowed counts with zero values

    Returns
    -------
    Tuple[Dict[bytes, int], Dict[bytes, int], List[bytes]]
        Tuple containing:
        - allowed_counts: Dictionary of counts for allowed UMIs
        - banned_counts: Dictionary of counts for UMIs not in allowed list
        - missing: List of UMIs in allowed list but not in counts

    Notes
    -----
    This function is commonly used in statistics calculations to filter data
    based on a predefined set of allowed UMIs and handle missing values.
    """
    allowed_counts = {k: v for k, v in counts.items() if k in allowed_list}
    banned_counts = {k: v for k, v in counts.items() if k not in allowed_list}
    missing = [umi for umi in allowed_list if umi not in counts]

    if add_missing:
        for umi in missing:
            allowed_counts[umi] = 0

    logger.debug(
        f"Split {len(counts)} counts: {len(allowed_counts)} allowed, {len(banned_counts)} banned, {len(missing)} missing"
    )
    return allowed_counts, banned_counts, missing


# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
