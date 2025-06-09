"""Single UMI library statistics"""

from typing import Dict, List, Optional, Union, Sequence
import numpy as np
import pandas as pd
from ..umi import UMI
from .base import BaseStatistic
from .utils import split_counts_by_allowed_list

class UMIStats(BaseStatistic):
    """Calculate UMI statistics
    
    A class with some helpful methods for calculating UMI statistics.
    """

    def __init__(self, umi: UMI, use_corrected: bool = True, allowed_list: Optional[List[str]] = None):
        """Initialize UMI stats calculator
        
        """
        super().__init__(umi)
        self.use_corrected = use_corrected
        self.allowed_list = allowed_list

    @property
    def _allowed_list(self) -> List[str]:
        """Get allowed list"""
        if self.allowed_list:
            return [umi.encode('ascii') for umi in self.allowed_list]

    @property
    def _counts(self) -> Dict[bytes, int]:
        """Get counts"""
        if self.use_corrected:
            return self.umi.corrected_counts
        return self.umi._counts 


class GiniCoefficient(UMIStats):
    """Calculate Gini coefficient for UMI counts"""
    
    def __init__(self, umi: UMI, use_corrected: bool = True, allowed_list: Optional[List[str]] = None):
        """Initialize Gini coefficient calculator
        
        Args:
            umi: UMI object to calculate statistics on
            use_corrected: If True, use corrected counts. If False, use original counts.
            allowed_list: Optional list of allowed keys. If provided, missing keys will be
                         treated as having zero counts in the calculation.
        """
        super().__init__(umi, use_corrected, allowed_list)
    
    @staticmethod
    def calculate_gini(counts: Sequence[float]) -> Optional[float]:
        """Calculate Gini coefficient from a sequence of counts
        
        Args:
            counts: Sequence of count values
            
        Returns:
            Gini coefficient or None if calculation is not possible
            
        Notes:
            The Gini coefficient is calculated using the formula:
            G = (2 * sum(i * y_i)) / (n * sum(y_i)) - (n + 1) / n
            where i is the rank and y_i is the count value
        """
        if not counts:
            return None
            
        # Sort counts in ascending order
        sorted_counts = sorted(counts)
        n = len(sorted_counts)
        
        # If all counts are zero, return None
        total = sum(sorted_counts)
        if total == 0:
            return None
            
        # Calculate the Lorenz curve
        index = list(range(1, n + 1))
        gini = ((2 * sum(i * y for i, y in zip(index, sorted_counts))) / 
                (n * total)) - ((n + 1) / n)
        
        return gini
    
    def run(self) -> Dict[str, Optional[float]]:
        """Calculate the Gini coefficient for the UMI object
        
        Returns:
            Dictionary containing the Gini coefficient
        """
        if self.allowed_list:
            counts, _, _ = split_counts_by_allowed_list(self._counts, self._allowed_list, add_missing=True)
        else:
            counts = self._counts

        return GiniCoefficient.calculate_gini(list(counts.values()))


class ShannonDiversity(UMIStats):
    """Calculate Shannon diversity index for UMI counts"""
    
    def __init__(self, umi: UMI, use_corrected: bool = True, base: float = 2.0, allowed_list: Optional[List[str]] = None):
        """Initialize Shannon diversity calculator
        
        Args:
            umi: UMI object to calculate statistics on
            use_corrected: If True, use corrected counts. If False, use original counts.
            base: Base of the logarithm (default: 2.0 for bits)
        """
        super().__init__(umi, use_corrected, allowed_list)
        self.base = base
    
    @staticmethod
    def calculate_shannon(counts: Sequence[float], base: float = 2.0) -> Optional[float]:
        """Calculate Shannon diversity index from a sequence of counts
        
        Args:
            counts: Sequence of count values
            base: Base of the logarithm (default: 2.0 for bits)
            
        Returns:
            Shannon diversity index or None if calculation is not possible
            
        Notes:
            The Shannon diversity index is calculated using the formula:
            H = -sum(p_i * log(p_i)) / log(base)
            where p_i is the proportion of each count
        """
        if not counts:
            return None
        
        # Convert to numpy array for efficient calculation
        counts_array = np.array([c for c in counts if c > 0])
        total = np.sum(counts_array)
        
        if total == 0:
            return None
        
        # Calculate proportions
        proportions = counts_array / total
        
        # Calculate Shannon diversity
        shannon = -np.sum(proportions * np.log(proportions) / np.log(base))
        
        return shannon
    
    def run(self) -> Dict[str, Optional[float]]:
        """Calculate the Shannon diversity index for the UMI object
        
        Returns:
            Dictionary containing the Shannon diversity index
        """
        
        if self.allowed_list:
            counts, _, _ = split_counts_by_allowed_list(self._counts, self._allowed_list, add_missing=False) 
        else:
            counts = self._counts
        
        return ShannonDiversity.calculate_shannon(list(counts.values()), self.base)


class SimpsonDiversity(UMIStats):
    """Calculate Simpson's diversity index for UMI counts"""
    
    def __init__(self, umi: UMI, use_corrected: bool = True, allowed_list: Optional[List[str]] = None):
        """Initialize Simpson's diversity calculator
        
        Args:
            umi: UMI object to calculate statistics on
            use_corrected: If True, use corrected counts. If False, use original counts.
            allowed_list: Optional list of allowed keys. If provided, missing keys will be
                         treated as having zero counts in the calculation.
        """
        super().__init__(umi, use_corrected, allowed_list)
    
    @staticmethod
    def calculate_simpson(counts: Sequence[float]) -> Optional[float]:
        """Calculate Simpson's diversity index from a sequence of counts
        
        Args:
            counts: Sequence of count values
            
        Returns:
            Simpson's diversity index or None if calculation is not possible
            
        Notes:
            Simpson's diversity index is calculated using the formula:
            D = 1 - sum(p_i^2)
            where p_i is the proportion of each count
        """
        if not counts:
            return None
        
        # Convert dict_values to list before numpy array conversion
        counts_array = np.array(list(counts))
        total = np.sum(counts_array)
        
        if total == 0:
            return None
        
        # Calculate proportions
        proportions = counts_array / total
        
        # Calculate Simpson's diversity (1 - D)
        simpson = 1 - np.sum(proportions ** 2)
        
        return simpson
    
    def run(self) -> Dict[str, Optional[float]]:
        """Calculate Simpson's diversity index for the UMI object
        
        Returns:
            Dictionary containing Simpson's diversity index
        """
        if self.allowed_list:
            counts, _, _ = split_counts_by_allowed_list(self._counts, self._allowed_list, add_missing=False) 
        else:
            counts = self._counts

        return SimpsonDiversity.calculate_simpson(list(counts.values()))


class UMIRecoveryRate(UMIStats):
    """Calculate UMI recovery rate"""
    
    def __init__(self, umi: UMI, use_corrected: bool = True, allowed_list: Optional[List[str]] = None):
        """Initialize UMI recovery rate calculator
        
        Args:
            umi: UMI object to calculate statistics on
            use_corrected: If True, use corrected counts. If False, use original counts.
            allowed_list: Optional list of allowed UMIs. If provided, recovery rate is
                         calculated as the ratio of observed allowed UMIs to total allowed UMIs.
                         If not provided, assumes exponential distribution and calculates
                         theoretical recovery rate.
        """
        super().__init__(umi, use_corrected, allowed_list)

    @staticmethod
    def calculate_recovery_rate_limited(counts: Dict[bytes, int], allowed_list: List[bytes]) -> Optional[float]:
        """Calculate UMI recovery rate from observed and total unique counts
        
        Args:
            counts: Dictionary of UMIs and their counts
            allowed_list: List of allowed UMIs
        """       
        _, _, missing = split_counts_by_allowed_list(counts, allowed_list, add_missing=True)
        return (len(allowed_list) - len(missing)) / len(allowed_list)
    
    def run(self) -> Dict[str, Optional[float]]:
        """Calculate UMI recovery rate for the UMI object
        
        Returns:
            Dictionary containing UMI recovery rate
        """
        if not self.allowed_list:
            return None
        
        return UMIRecoveryRate.calculate_recovery_rate_limited(self._counts, self._allowed_list)


class UMIEfficiencyRate(UMIStats):
    """Calculate UMI efficiency rate (fraction of reads contributing to allowed UMIs)"""
    
    def __init__(self, umi: UMI, use_corrected: bool = True, allowed_list: Optional[List[str]] = None):
        """Initialize UMI efficiency rate calculator
        
        Args:
            umi: UMI object to calculate statistics on
            allowed_list: List of allowed UMIs. Efficiency rate is calculated as the
                         fraction of reads that contributed to these allowed UMIs.
            use_corrected: If True, use corrected counts. If False, use original counts.
        """ 
        super().__init__(umi, use_corrected, allowed_list)
    
    @staticmethod
    def calculate_efficiency_rate(counts: Dict[bytes, int], allowed_list: List[bytes]) -> Optional[float]:
        """Calculate UMI efficiency rate from allowed and total reads
        
        Args:
            counts: Dictionary of UMIs and their counts
            allowed_list: List of allowed UMIs
            
        Returns:
            UMI efficiency rate or None if calculation is not possible
            
        Notes:
            UMI efficiency rate is calculated as allowed_reads / total_reads
        """
        
        if not counts:
            return None
        
        counts, banned, _ = split_counts_by_allowed_list(counts, allowed_list, add_missing=True)
        
        wanted_reads = sum(counts.values())
        total_reads = sum(counts.values()) + sum(banned.values())

        if total_reads == 0:
            return None

        return wanted_reads / total_reads
    
    def run(self) -> Dict[str, Optional[float]]:
        """Calculate UMI efficiency rate for the UMI object
        
        Returns:
            Dictionary containing UMI efficiency rate
        """
        if not self.allowed_list:
            return None
        
        return UMIEfficiencyRate.calculate_efficiency_rate(self._counts, self._allowed_list)
        

class UMIErrorRate(BaseStatistic):
    """Calculate UMI error rate based on mismatches between original and corrected UMIs"""
    
    def __init__(self, umi: UMI):
        """Initialize UMI error rate calculator
        
        Args:
            umi: UMI object to calculate statistics on
        """
        super().__init__(umi)
    
    @staticmethod
    def hamming_distance(seq1: bytes, seq2: bytes) -> int:
        """Calculate Hamming distance between two sequences
        
        Args:
            seq1: First sequence
            seq2: Second sequence
            
        Returns:
            Number of mismatches between sequences
        """
        if len(seq1) != len(seq2):
            raise ValueError("Sequences must be of the same length")
        return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))
    
    @staticmethod
    def calculate_error_rate(mapping: Dict[bytes, bytes], counts: Dict[bytes, int]) -> Optional[float]:
        """Calculate UMI error rate from mismatches and total reads
        
        Args:
            mapping: Dictionary of original and corrected UMIs
            counts: Dictionary of UMIs and their counts
            
        Returns:
            UMI error rate or None if calculation is not possible
            
        Notes:
            UMI error rate is calculated as total mismatches / total reads
        """
        # Calculate total mismatches
        if not mapping or not counts:
            return None
        
        total_mismatches = 0
        for original, corrected in mapping.items():
            if original != corrected:  # Only calculate mismatches if UMI was corrected
                mismatches = UMIErrorRate.hamming_distance(original, corrected)
                total_mismatches += mismatches * counts[original]
        
        total_reads = sum(counts.values())

        return total_mismatches / total_reads
    
    def run(self) -> Dict[str, Optional[float]]:
        """Calculate UMI error rate for the UMI object
        
        Returns:
            Dictionary containing UMI error rate
        """
        if not self.umi._mapping or not self.umi._counts:
            return None
        
        return UMIErrorRate.calculate_error_rate(self.umi._mapping, self.umi._counts)


class UMIRedundancy(UMIStats):
    """Calculate UMI redundancy (average reads per unique UMI)"""
    
    def __init__(self, umi: UMI, use_corrected: bool = True, allowed_list: Optional[List[str]] = None):
        """Initialize UMI redundancy calculator
        
        Args:
            umi: UMI object to calculate statistics on
            use_corrected: If True, use corrected counts. If False, use original counts.
        """
        super().__init__(umi, use_corrected, allowed_list)
    
    @staticmethod
    def calculate_redundancy(counts: Dict[bytes, int]) -> Optional[float]:
        """Calculate UMI redundancy from total and unique counts
        
        Args:
            counts: Dictionary of UMIs and their counts
            
        Returns:
            UMI redundancy or None if calculation is not possible
            
        Notes:
            UMI redundancy is calculated as total_reads / unique_umis
        """
        total_reads = sum(counts.values())
        unique_umis = len(counts)

        if unique_umis == 0:
            return None
        
        return total_reads / unique_umis
    
    def run(self) -> Dict[str, Optional[float]]:
        """Calculate UMI redundancy for the UMI object
        
        Returns:
            Dictionary containing UMI redundancy
        """
        
        if self.allowed_list:
            counts, _, _ = split_counts_by_allowed_list(self._counts, self._allowed_list, add_missing=False) 
        else:
            counts = self._counts
        
        if not counts:
            return None
        
        return UMIRedundancy.calculate_redundancy(counts) 
