"""Pairwise and UMI-level comparison statistics"""

from typing import Dict, List, Optional, Union, Sequence
import numpy as np
import pandas as pd
from scipy import stats
from ..umi import UMI
from .base import BasePairwiseStatistic, BaseDifferentialStatistic
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
            Dictionary containing Jaccard similarity coefficient
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
            Dictionary containing Bray-Curtis dissimilarity coefficient
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
            Dictionary containing Spearman correlation coefficient
        """
        counts1 = self.umi1.corrected_counts if self.use_corrected else self.umi1._counts
        counts2 = self.umi2.corrected_counts if self.use_corrected else self.umi2._counts
        
        if self.allowed_list:
            counts1, _, _ = split_counts_by_allowed_list(counts1, self._allowed_list, add_missing=False)
            counts2, _, _ = split_counts_by_allowed_list(counts2, self._allowed_list, add_missing=False)
        
        return SpearmanCorrelation.calculate_spearman(counts1, counts2)


class DifferentialAbundance(BaseDifferentialStatistic):
    """Calculate differential abundance scores between UMI groups"""
    
    def __init__(self, umis: List[UMI], groups: List[str], use_corrected: bool = True, allowed_list: Optional[List[str]] = None):
        """Initialize differential abundance calculator
        
        Args:
            umis: List of UMI objects
            groups: List of group labels corresponding to UMI objects
            use_corrected: If True, use corrected counts. If False, use original counts.
            allowed_list: Optional list of allowed UMIs. If provided, only these UMIs will be considered.
        """
        super().__init__(umis, groups)
        self.use_corrected = use_corrected
        self.allowed_list = allowed_list

    @property
    def _allowed_list(self) -> Optional[List[bytes]]:
        """Get allowed list in bytes format"""
        if self.allowed_list:
            return [umi.encode('ascii') for umi in self.allowed_list]
        return None
    
    @staticmethod
    def calculate_differential_abundance(counts_df: pd.DataFrame, group_map: Dict[str, str]) -> pd.DataFrame:
        """Calculate differential abundance scores for each UMI
        
        Args:
            counts_df: DataFrame where rows are UMIs and columns are samples
            group_map: Dictionary mapping sample names to group labels
            
        Returns:
            DataFrame with differential abundance metrics for each UMI
            
        Raises:
            ValueError: If group_map does not contain exactly 2 groups
        """

        # Validate counts_df
        if len(counts_df.index) == 0:
            raise ValueError("Counts DataFrame is empty")

        # Validate group map
        unique_groups = set(group_map.values())
        if len(unique_groups) != 2:
            raise ValueError("Differential abundance calculation requires exactly 2 groups")
            
        # Get group names
        group1, group2 = sorted(unique_groups)
        
        # Get samples for each group
        group1_samples = [s for s, g in group_map.items() if g == group1]
        group2_samples = [s for s, g in group_map.items() if g == group2]
        
        # Calculate metrics for each UMI
        results = []
        for umi in counts_df.index:
            # Get counts for each group
            group1_counts = counts_df.loc[umi, group1_samples].values
            group2_counts = counts_df.loc[umi, group2_samples].values
            
            # Calculate log2 fold change
            mean1 = np.mean(group1_counts)
            mean2 = np.mean(group2_counts)
            
            if mean1 == 0:
                lfc = float('inf') if mean2 > 0 else 0.0
            else:
                lfc = np.log2(mean2 / mean1)
            
            # Calculate effect size (Cohen's d)
            std1 = np.std(group1_counts)
            std2 = np.std(group2_counts)
            denominator = np.sqrt((std1**2 + std2**2) / 2)
            if denominator == 0:
                effect_size = 0.0
            else:
                effect_size = (mean2 - mean1) / denominator
            
            # Calculate p-value using Mann-Whitney U test
            _, p_value = stats.mannwhitneyu(group1_counts, group2_counts, alternative='two-sided')
            
            results.append({
                'umi': umi,
                'log2_fold_change': lfc,
                'effect_size': effect_size,
                'p_value': p_value,
                'mean_group1': mean1,
                'mean_group2': mean2,
                'std_group1': std1,
                'std_group2': std2
            })
            #print(umi, lfc, effect_size, p_value)
            
        return pd.DataFrame(results).set_index('umi')
    
    def run(self) -> pd.DataFrame:
        """Calculate differential abundance scores between UMI groups
        
        Returns:
            DataFrame containing differential abundance metrics for each UMI
        """
        # Create DataFrame of counts
        all_umis = set()
        for umi in self.umis:
            counts = umi.corrected_counts if self.use_corrected else umi._counts
            if self.allowed_list:
                counts, _, _ = split_counts_by_allowed_list(counts, self._allowed_list, add_missing=True)
            all_umis.update(counts.keys())
            
        # Create DataFrame with all UMIs and samples
        counts_df = pd.DataFrame(index=list(all_umis))
        for i, (umi, group) in enumerate(zip(self.umis, self.groups)):
            counts = umi.corrected_counts if self.use_corrected else umi._counts
            if self.allowed_list:
                counts, _, _ = split_counts_by_allowed_list(counts, self._allowed_list, add_missing=True)
            counts_df[f'sample_{i}'] = pd.Series(counts)
            
        # Create group mapping
        group_map = {f'sample_{i}': group for i, group in enumerate(self.groups)}
            
        return DifferentialAbundance.calculate_differential_abundance(counts_df, group_map) 