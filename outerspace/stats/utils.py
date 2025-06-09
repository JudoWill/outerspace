"""Utility functions for stats"""

from typing import Dict, List, Tuple

def split_counts_by_allowed_list(counts: Dict, allowed_list: List, add_missing: bool = True) -> Tuple[Dict, Dict, List]:
    """Split counts by allowed list
    
    Args:
        counts: Dictionary of counts
        allowed_list: List of allowed UMIs
        add_missing: Whether to add missing UMIs to the counts
        
    Returns
        - allowed_counts: Dictionary of counts for allowed UMIs
        - banned_counts: Dictionary of counts for UMIs not in allowed list
        - missing : List of UMIs in allowed list but not in counts
    """
    allowed_counts = {k: v for k, v in counts.items() if k in allowed_list}
    banned_counts = {k: v for k, v in counts.items() if k not in allowed_list}
    missing = [umi for umi in allowed_list if umi not in counts]
    if add_missing:
        for umi in missing:
            allowed_counts[umi] = 0
    return allowed_counts, banned_counts, missing