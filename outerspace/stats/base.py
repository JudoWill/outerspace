"""Base classes for UMI statistics calculations"""

from typing import Dict, List, Optional, TypeVar, Union
from abc import ABC, abstractmethod
from ..umi import UMI, UmiCollection

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
    def run(self) -> T:
        """Run the statistical calculation
        
        Returns:
            Calculated statistic value
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
    def calculate(cls, umi: UMI, **kwargs) -> T:
        """Create an instance and run the calculation
        
        Args:
            umi: UMI object to calculate statistics on
            **kwargs: Additional arguments specific to the statistic
            
        Returns:
            Calculated statistic value
        """
        instance = cls(umi, **kwargs)
        instance.validate()
        instance.normalize()
        return instance.run()

    @classmethod
    def calculate_collection(cls, collection: UmiCollection, **kwargs) -> Dict[str, T]:
        """Calculate statistics for each UMI in a collection
        
        Args:
            collection: UmiCollection object containing multiple UMIs
            **kwargs: Additional arguments specific to the statistic
            
        Returns:
            Dictionary mapping sample names to calculated statistics
        """
        results = {}
        for sample_name, umi in collection.umis.items():
            results[sample_name] = cls.calculate(umi, **kwargs)
        return results


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
    def run(self) -> T:
        """Run the pairwise statistical calculation
        
        Returns:
            Calculated statistic value
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
    def calculate(cls, umi1: UMI, umi2: UMI, **kwargs) -> T:
        """Create an instance and run the calculation
        
        Args:
            umi1: First UMI object
            umi2: Second UMI object
            **kwargs: Additional arguments specific to the statistic
            
        Returns:
            Calculated statistic value
        """
        instance = cls(umi1, umi2, **kwargs)
        instance.validate()
        instance.normalize()
        return instance.run()

    @classmethod
    def calculate_collection(cls, collection: UmiCollection, **kwargs) -> Dict[tuple[str, str], T]:
        """Calculate pairwise statistics for all pairs of UMIs in a collection
        
        Args:
            collection: UmiCollection object containing multiple UMIs
            **kwargs: Additional arguments specific to the statistic
            
        Returns:
            Dictionary mapping (sample1, sample2) tuples to calculated statistics
        """
        results = {}
        samples = list(collection.umis.keys())
        for i, sample1 in enumerate(samples):
            for sample2 in samples[i+1:]:
                results[(sample1, sample2)] = cls.calculate(
                    collection.umis[sample1],
                    collection.umis[sample2],
                    **kwargs
                )
        return results


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
    def run(self) -> T:
        """Run the differential statistical calculation
        
        Returns:
            Calculated statistic value
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
    def calculate(cls, umis: List[UMI], groups: List[str], **kwargs) -> T:
        """Create an instance and run the calculation
        
        Args:
            umis: List of UMI objects
            groups: List of group labels corresponding to UMI objects
            **kwargs: Additional arguments specific to the statistic
            
        Returns:
            Calculated statistic value
        """
        instance = cls(umis, groups, **kwargs)
        instance.validate()
        instance.normalize()
        return instance.run()

    @classmethod
    def calculate_collection(cls, collection: UmiCollection, groups: Dict[str, str], **kwargs) -> T:
        """Calculate differential statistics for a collection of UMIs
        
        Args:
            collection: UmiCollection object containing multiple UMIs
            groups: Dictionary mapping sample names to group labels
            **kwargs: Additional arguments specific to the statistic
            
        Returns:
            Calculated statistic value
            
        Raises:
            ValueError: If any sample in collection is not in groups
        """
        # Validate that all samples have group assignments
        missing = set(collection.umis.keys()) - set(groups.keys())
        if missing:
            raise ValueError(f"Samples missing group assignments: {missing}")
            
        # Create ordered lists of UMIs and groups
        umis = []
        group_labels = []
        for sample in collection.umis.keys():
            umis.append(collection.umis[sample])
            group_labels.append(groups[sample])
            
        return cls.calculate(umis, group_labels, **kwargs) 