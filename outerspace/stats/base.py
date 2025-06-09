"""Base classes for UMI statistics calculations"""

from typing import Dict, List, Optional, TypeVar
from abc import ABC, abstractmethod
from ..umi import UMI

T = TypeVar('T')

class BaseStatistic(ABC):
    """Base class for all statistics calculations"""
    
    def __init__(self, umi: UMI):
        """Initialize the statistic calculator
        
        Args:
            umi: UMI object to calculate statistics on
            **kwargs: Additional arguments specific to the statistic
        """
        self.umi = umi
    
    @abstractmethod
    def run(self) -> Dict[str, T]:
        """Run the statistical calculation
        
        Returns:
            Dictionary of calculated statistics
        """
        pass
    
    def validate(self) -> None:
        """Validate input parameters
        
        Raises:
            ValueError: If validation fails
        """
        pass
    
    def normalize(self) -> None:
        """Normalize data before calculation if needed"""
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
        instance.validate()
        instance.normalize()
        return instance.run()


class BasePairwiseStatistic(ABC):
    """Base class for pairwise statistics"""
    
    def __init__(self, umi1: UMI, umi2: UMI):
        """Initialize the pairwise statistic calculator
        
        Args:
            umi1: First UMI object
            umi2: Second UMI object
            **kwargs: Additional arguments specific to the statistic
        """
        self.umi1 = umi1
        self.umi2 = umi2
        self._result: Optional[Dict[str, T]] = None
    
    @abstractmethod
    def run(self) -> Dict[str, T]:
        """Run the pairwise statistical calculation
        
        Returns:
            Dictionary of calculated statistics
        """
        pass
    
    def validate(self) -> None:
        """Validate input parameters
        
        Raises:
            ValueError: If validation fails
        """
        pass
    
    def normalize(self) -> None:
        """Normalize data before calculation if needed"""
        pass
    
    @classmethod
    def calculate(cls, umi1: UMI, umi2: UMI, **kwargs) -> Dict[str, T]:
        """Create an instance and run the calculation
        
        Args:
            umi1: First UMI object
            umi2: Second UMI object
            **kwargs: Additional arguments specific to the statistic
            
        Returns:
            Dictionary of calculated statistics
        """
        instance = cls(umi1, umi2, **kwargs)
        instance.validate()
        instance.normalize()
        return instance.run()


class BaseDifferentialStatistic(ABC):
    """Base class for differential analysis statistics"""
    
    def __init__(self, umis: List[UMI], groups: List[str]):
        """Initialize the differential statistic calculator
        
        Args:
            umis: List of UMI objects
            groups: List of group labels corresponding to UMI objects
            **kwargs: Additional arguments specific to the statistic
        """
        if len(umis) != len(groups):
            raise ValueError("Number of UMIs must match number of groups")
        self.umis = umis
        self.groups = groups
        self._result: Optional[Dict[str, T]] = None
    
    @abstractmethod
    def run(self) -> Dict[str, T]:
        """Run the differential statistical calculation
        
        Returns:
            Dictionary of calculated statistics
        """
        pass
    
    def validate(self) -> None:
        """Validate input parameters
        
        Raises:
            ValueError: If validation fails
        """
        pass
    
    def normalize(self) -> None:
        """Normalize data before calculation if needed"""
        pass
    
    @classmethod
    def calculate(cls, umis: List[UMI], groups: List[str], **kwargs) -> Dict[str, T]:
        """Create an instance and run the calculation
        
        Args:
            umis: List of UMI objects
            groups: List of group labels corresponding to UMI objects
            **kwargs: Additional arguments specific to the statistic
            
        Returns:
            Dictionary of calculated statistics
        """
        instance = cls(umis, groups, **kwargs)
        instance.validate()
        instance.normalize()
        return instance.run() 