"""Pairwise and UMI-level comparison statistics"""

from typing import Dict, List, Optional, Union, Sequence
import numpy as np
import pandas as pd
from ..umi import UMI, UmiCollection
from .base import BasePairwiseStatistic
from .utils import split_counts_by_allowed_list


class JaccardSimilarity(BasePairwiseStatistic):
    """Calculate Jaccard similarity between two UMI libraries"""
    
    def __init__(self, umi1: UMI, umi2: UMI, use_corrected: bool = True, allowed_list: Optional[List[str]] = None):
        """Initialize Jaccard similarity calculator
        
        Args:
            umi1: First UMI object
            umi2: Second UMI object
            use_corrected: If True, use corrected counts. If False, use original counts.
            allowed_list: Optional list of allowed UMIs. If provided, only these UMIs will be considered.
        """
        super().__init__(umi1, umi2)
        self.use_corrected = use_corrected
        self.allowed_list = allowed_list

    @property
    def _allowed_list(self) -> Optional[List[bytes]]:
        """Get allowed list in bytes format"""
        if self.allowed_list:
            return [umi.encode('ascii') for umi in self.allowed_list]
        return None
    
    @staticmethod
    def calculate_jaccard(set1: set, set2: set) -> float:
        """Calculate Jaccard similarity between two sets
        
        Args:
            set1: First set
            set2: Second set
            
        Returns:
            Jaccard similarity coefficient
        """
        intersection = len(set1.intersection(set2))
        union = len(set1.union(set2))
        return intersection / union if union > 0 else 0.0
    
    def run(self) -> float:
        """Calculate Jaccard similarity between UMI libraries
        
        Returns:
            Jaccard similarity coefficient
        """
        counts1 = self.umi1.corrected_counts if self.use_corrected else self.umi1._counts
        counts2 = self.umi2.corrected_counts if self.use_corrected else self.umi2._counts
        
        if self.allowed_list:
            counts1, _, _ = split_counts_by_allowed_list(counts1, self._allowed_list, add_missing=False)
            counts2, _, _ = split_counts_by_allowed_list(counts2, self._allowed_list, add_missing=False)
        
        set1 = set(counts1.keys())
        set2 = set(counts2.keys())
        
        return JaccardSimilarity.calculate_jaccard(set1, set2)


class BrayCurtisDissimilarity(BasePairwiseStatistic):
    """Calculate Bray-Curtis dissimilarity between two UMI libraries"""
    
    def __init__(self, umi1: UMI, umi2: UMI, use_corrected: bool = True, allowed_list: Optional[List[str]] = None):
        """Initialize Bray-Curtis dissimilarity calculator
        
        Args:
            umi1: First UMI object
            umi2: Second UMI object
            use_corrected: If True, use corrected counts. If False, use original counts.
            allowed_list: Optional list of allowed UMIs. If provided, only these UMIs will be considered.
        """
        super().__init__(umi1, umi2)
        self.use_corrected = use_corrected
        self.allowed_list = allowed_list

    @property
    def _allowed_list(self) -> Optional[List[bytes]]:
        """Get allowed list in bytes format"""
        if self.allowed_list:
            return [umi.encode('ascii') for umi in self.allowed_list]
        return None
    
    @staticmethod
    def calculate_bray_curtis(counts1: Dict[bytes, int], counts2: Dict[bytes, int]) -> float:
        """Calculate Bray-Curtis dissimilarity between two count distributions
        
        Args:
            counts1: First count distribution
            counts2: Second count distribution
            
        Returns:
            Bray-Curtis dissimilarity coefficient
        """
        # Get all unique keys
        all_keys = set(counts1.keys()) | set(counts2.keys())
        
        # Calculate sums
        sum1 = sum(counts1.values())
        sum2 = sum(counts2.values())
        
        if sum1 == 0 and sum2 == 0:
            return 0.0
            
        # Calculate differences and sums
        diff_sum = 0
        total_sum = 0
        
        for key in all_keys:
            count1 = counts1.get(key, 0)
            count2 = counts2.get(key, 0)
            diff_sum += abs(count1 - count2)
            total_sum += count1 + count2
            
        return diff_sum / total_sum if total_sum > 0 else 0.0
    
    def run(self) -> float:
        """Calculate Bray-Curtis dissimilarity between UMI libraries
        
        Returns:
            Bray-Curtis dissimilarity coefficient
        """
        counts1 = self.umi1.corrected_counts if self.use_corrected else self.umi1._counts
        counts2 = self.umi2.corrected_counts if self.use_corrected else self.umi2._counts
        
        if self.allowed_list:
            counts1, _, _ = split_counts_by_allowed_list(counts1, self._allowed_list, add_missing=False)
            counts2, _, _ = split_counts_by_allowed_list(counts2, self._allowed_list, add_missing=False)
        
        return BrayCurtisDissimilarity.calculate_bray_curtis(counts1, counts2)


class FoldChange(BasePairwiseStatistic):
    """Calculate fold change between two UMI libraries"""
    
    def __init__(self, umi1: UMI, umi2: UMI, use_corrected: bool = True, allowed_list: Optional[List[str]] = None):
        """Initialize fold change calculator
        
        Args:
            umi1: First UMI object
            umi2: Second UMI object
            use_corrected: If True, use corrected counts. If False, use original counts.
            allowed_list: Optional list of allowed UMIs. If provided, only these UMIs will be considered.
        """
        super().__init__(umi1, umi2)
        self.use_corrected = use_corrected
        self.allowed_list = allowed_list

    @property
    def _allowed_list(self) -> Optional[List[bytes]]:
        """Get allowed list in bytes format"""
        if self.allowed_list:
            return [umi.encode('ascii') for umi in self.allowed_list]
        return None
    
    @staticmethod
    def calculate_fold_change(counts1: Dict[bytes, int], counts2: Dict[bytes, int]) -> Dict[bytes, float]:
        """Calculate fold change for each UMI between two count distributions
        
        Args:
            counts1: First count distribution
            counts2: Second count distribution
            
        Returns:
            Dictionary mapping UMIs to their fold changes
        """
        fold_changes = {}
        all_keys = set(counts1.keys()) | set(counts2.keys())
        
        for key in all_keys:
            count1 = counts1.get(key, 0)
            count2 = counts2.get(key, 0)
            
            # Avoid division by zero
            if count1 == 0:
                fold_changes[key] = float('inf') if count2 > 0 else 0.0
            else:
                fold_changes[key] = count2 / count1
                
        return fold_changes
    
    def run(self) -> Dict[bytes, float]:
        """Calculate fold changes between UMI libraries
        
        Returns:
            Dictionary containing fold changes for each UMI
        """
        counts1 = self.umi1.corrected_counts if self.use_corrected else self.umi1._counts
        counts2 = self.umi2.corrected_counts if self.use_corrected else self.umi2._counts
        
        if self.allowed_list:
            counts1, _, _ = split_counts_by_allowed_list(counts1, self._allowed_list, add_missing=False)
            counts2, _, _ = split_counts_by_allowed_list(counts2, self._allowed_list, add_missing=False)
        
        return FoldChange.calculate_fold_change(counts1, counts2)


class SpearmanCorrelation(BasePairwiseStatistic):
    """Calculate Spearman correlation between two UMI libraries"""
    
    def __init__(self, umi1: UMI, umi2: UMI, use_corrected: bool = True, allowed_list: Optional[List[str]] = None):
        """Initialize Spearman correlation calculator
        
        Args:
            umi1: First UMI object
            umi2: Second UMI object
            use_corrected: If True, use corrected counts. If False, use original counts.
            allowed_list: Optional list of allowed UMIs. If provided, only these UMIs will be considered.
        """
        super().__init__(umi1, umi2)
        self.use_corrected = use_corrected
        self.allowed_list = allowed_list

    @property
    def _allowed_list(self) -> Optional[List[bytes]]:
        """Get allowed list in bytes format"""
        if self.allowed_list:
            return [umi.encode('ascii') for umi in self.allowed_list]
        return None
    
    @staticmethod
    def calculate_spearman(counts1: Dict[bytes, int], counts2: Dict[bytes, int]) -> float:
        """Calculate Spearman correlation between two count distributions
        
        Args:
            counts1: First count distribution
            counts2: Second count distribution
            
        Returns:
            Spearman correlation coefficient
        """
        # Get common keys
        common_keys = set(counts1.keys()) & set(counts2.keys())
        
        if not common_keys:
            return 0.0
            
        # Create arrays of counts for common keys
        counts1_array = np.array([counts1[key] for key in common_keys])
        counts2_array = np.array([counts2[key] for key in common_keys])
        
        # Calculate Spearman correlation
        correlation = pd.Series(counts1_array).corr(pd.Series(counts2_array), method='spearman')
        
        return correlation
    
    def run(self) -> float:
        """Calculate Spearman correlation between UMI libraries
        
        Returns:
            Spearman correlation coefficient
        """
        counts1 = self.umi1.corrected_counts if self.use_corrected else self.umi1._counts
        counts2 = self.umi2.corrected_counts if self.use_corrected else self.umi2._counts
        
        if self.allowed_list:
            counts1, _, _ = split_counts_by_allowed_list(counts1, self._allowed_list, add_missing=False)
            counts2, _, _ = split_counts_by_allowed_list(counts2, self._allowed_list, add_missing=False)
        
        return SpearmanCorrelation.calculate_spearman(counts1, counts2)


