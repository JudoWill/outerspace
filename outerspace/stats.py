"""Module for statistical calculations on UMI data"""

from typing import Dict, List, Optional, Union, TypeVar, Generic
from abc import ABC, abstractmethod
from .umi import UMI

T = TypeVar('T')

class BaseStatistic(ABC, Generic[T]):
    """Base class for all statistics calculations"""
    
    def __init__(self, umi: UMI, **kwargs):
        """Initialize the statistic calculator
        
        Args:
            umi: UMI object to calculate statistics on
            **kwargs: Additional arguments specific to the statistic
        """
        self.umi = umi
        self.kwargs = kwargs
    
    @abstractmethod
    def run(self) -> Dict[str, T]:
        """Run the statistical calculation
        
        Returns:
            Dictionary of calculated statistics
        """
        pass
    
    @classmethod
    def calculate(cls, umi: UMI, **kwargs) -> Dict[str, T]:
        """Create an instance and run the calculation
        
        Args:
            umi: UMI object to calculate statistics on
            **kwargs: Additional arguments specific to the statistic
            
        Returns:
            Dictionary of calculated statistics
        """
        instance = cls(umi, **kwargs)
        return instance.run()


class GiniCoefficient(BaseStatistic[Optional[float]]):
    """Calculate Gini coefficient for UMI counts"""
    
    def __init__(self, umi: UMI, use_corrected: bool = True, allowed_list: Optional[List[str]] = None):
        """Initialize Gini coefficient calculator
        
        Args:
            umi: UMI object to calculate statistics on
            use_corrected: If True, use corrected counts. If False, use original counts.
            allowed_list: Optional list of allowed keys. If provided, missing keys will be
                         treated as having zero counts in the calculation.
        """
        super().__init__(umi, use_corrected=use_corrected, allowed_list=allowed_list)
    
    def run(self) -> Dict[str, Optional[float]]:
        """Calculate the Gini coefficient
        
        Returns:
            Dictionary containing the Gini coefficient
        """
        counts = self.umi.corrected_counts if self.kwargs['use_corrected'] else self.umi._counts
        
        if self.kwargs['allowed_list']:
            # Create a dictionary with all allowed keys, using 0 for missing ones
            # Convert string keys to bytes for lookup
            all_counts = {key.encode('ascii'): counts.get(key.encode('ascii'), 0) 
                         for key in self.kwargs['allowed_list']}
        else:
            all_counts = counts
            
        if not all_counts:
            return {'gini_coefficient': None}
        
        # Sort counts in ascending order
        sorted_counts = sorted(all_counts.values())
        n = len(sorted_counts)
        
        # If all counts are zero, return None (no data to measure inequality)
        total = sum(sorted_counts)
        if total == 0:
            return {'gini_coefficient': None}
            
        # Calculate the Lorenz curve
        index = list(range(1, n + 1))
        gini = ((2 * sum(i * y for i, y in zip(index, sorted_counts))) / 
                (n * total)) - ((n + 1) / n)
        
        return {'gini_coefficient': gini}
