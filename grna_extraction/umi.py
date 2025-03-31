"""Module for handling UMI clustering and correction"""

from typing import Dict, List, Optional, Union
from collections import Counter
from umi_tools import UMIClusterer
import csv
import os
from pathlib import Path


class UMI:
    """Class for handling UMI clustering and correction"""
    
    def __init__(self, mismatches: int = 2, method: str = "adjacency"):
        """Initialize UMI clusterer
        
        Args:
            mismatches: Number of mismatches allowed for clustering
            method: Clustering method to use (cluster, adjacency, or directional)
        """
        self.clusterer = UMIClusterer(cluster_method=method)
        self.mismatches = mismatches
        self._counts: Dict[bytes, int] = {}
        self._mapping: Dict[bytes, bytes] = {}
        self._corrected_counts: Optional[Dict[bytes, int]] = None
    
    def consume(self, umi: Union[str, bytes]) -> None:
        """Add a UMI to the counts dictionary
        
        Args:
            umi: UMI sequence to add
        """
        if isinstance(umi, str):
            umi = umi.encode('ascii')
        self._counts[umi] = self._counts.get(umi, 0) + 1
        self._corrected_counts = None  # Reset corrected counts
    
    def create_mapping(self) -> None:
        """Create mapping between original and corrected barcodes"""
        if not self._counts:
            return
            
        # Create clusters
        clusters = self.clusterer(self._counts, self.mismatches)
        
        # Create mapping
        self._mapping = {}
        for cluster in clusters:
            key = cluster[0]  # Use first barcode in cluster as representative
            for item in cluster:
                self._mapping[item] = key
    
    def __getitem__(self, umi: Union[str, bytes]) -> bytes:
        """Get corrected UMI for a given UMI
        
        Args:
            umi: Original UMI sequence
            
        Returns:
            Corrected UMI sequence
        """
        if isinstance(umi, str):
            umi = umi.encode('ascii')
        return self._mapping.get(umi, umi)
    
    def get(self, umi: Union[str, bytes], default: Optional[bytes] = None) -> Optional[bytes]:
        """Get corrected UMI for a given UMI with optional default
        
        Args:
            umi: Original UMI sequence
            default: Default value if UMI not found
            
        Returns:
            Corrected UMI sequence or default value
        """
        if isinstance(umi, str):
            umi = umi.encode('ascii')
        return self._mapping.get(umi, default)
    
    @property
    def corrected_counts(self) -> Dict[bytes, int]:
        """Get counts of corrected barcodes"""
        if self._corrected_counts is None:
            self._corrected_counts = {}
            for orig_bc, count in self._counts.items():
                corrected = self._mapping.get(orig_bc, orig_bc)
                self._corrected_counts[corrected] = self._corrected_counts.get(corrected, 0) + count
        return self._corrected_counts
    
    @classmethod
    def from_csv(cls, filepath: Union[str, Path], column: str, 
                 mismatches: int = 2, method: str = "adjacency",
                 sep: str = ",", correct: bool = True) -> 'UMI':
        """Create UMI object from CSV file
        
        Args:
            filepath: Path to CSV file
            column: Column containing UMIs
            mismatches: Number of mismatches allowed for clustering
            method: Clustering method to use
            sep: CSV separator
            correct: If True, perform clustering. If False, assume data is pre-clustered.
            
        Returns:
            UMI object
        """
        umi = cls(mismatches=mismatches, method=method)
        
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f, delimiter=sep)
            if column not in reader.fieldnames:
                raise ValueError(f"Column {column} not found in CSV file")
            
            for row in reader:
                umi.consume(row[column])
        
        if correct:
            umi.create_mapping()
        else:
            # For pre-clustered data, use each UMI as its own cluster
            umi._mapping = {bc: bc for bc in umi._counts.keys()}
            umi._corrected_counts = umi._counts.copy()
        
        return umi
    
    @classmethod
    def from_csvs(cls, directory: Union[str, Path], column: str,
                  mismatches: int = 2, method: str = "adjacency",
                  sep: str = ",", correct: bool = True) -> 'UMI':
        """Create UMI object from multiple CSV files
        
        Args:
            directory: Directory containing CSV files
            column: Column containing UMIs
            mismatches: Number of mismatches allowed for clustering
            method: Clustering method to use
            sep: CSV separator
            correct: If True, perform clustering. If False, assume data is pre-clustered.
            
        Returns:
            UMI object
        """
        umi = cls(mismatches=mismatches, method=method)
        directory = Path(directory)
        
        for csv_file in directory.glob("*.csv"):
            with open(csv_file, 'r') as f:
                reader = csv.DictReader(f, delimiter=sep)
                if column not in reader.fieldnames:
                    print(f"Warning: Column {column} not found in {csv_file}")
                    continue
                
                for row in reader:
                    umi.consume(row[column])
        
        if correct:
            umi.create_mapping()
        else:
            # For pre-clustered data, use each UMI as its own cluster
            umi._mapping = {bc: bc for bc in umi._counts.keys()}
            umi._corrected_counts = umi._counts.copy()
        
        return umi 