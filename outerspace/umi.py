"""Module for handling UMI clustering and correction"""

from typing import Dict, List, Optional, Union, Tuple
from collections import Counter
from umi_tools import UMIClusterer
import csv
import os
from pathlib import Path
import sys
import logging


def _cluster_umis(counts: Dict[bytes, int], 
                 mismatches: int = 2, 
                 method: str = "adjacency",
                 allowed_list: Optional[List[bytes]] = None) -> Tuple[Dict[bytes, int], Dict[bytes, bytes]]:
    """Cluster UMIs and return collapsed counts and mapping
    
    Args:
        counts: Dictionary mapping UMIs to their counts
        mismatches: Number of mismatches allowed for clustering
        method: Clustering method to use
        allowed_list: Optional list of allowed UMIs. If provided, highest priority UMI in cluster will be used as representative
        
    Returns:
        Tuple of (collapsed_counts, mapping)
    """
    if not counts:
        return {}, {}
        
    # If no mismatches allowed, each barcode maps to itself
    if mismatches == 0:
        return counts.copy(), {bc: bc for bc in counts.keys()}
        
    # Create clusters
    clusterer = UMIClusterer(cluster_method=method)
    clusters = clusterer(counts, mismatches)
    
    # Create mapping
    mapping = {}
    collapsed_counts = {}
    
    for cluster in clusters:
        if allowed_list:
            # Find highest priority UMI in cluster that's in allowed_list
            key = None
            # Find the UMI in the cluster that is in the allowed list
            # and has the highest count
            for bc in cluster:
                if bc in allowed_list:
                    if key is None or counts[bc] > counts[key]:
                        key = bc
            if key is None:
                # If no allowed UMIs in cluster, use most common
                key = max(cluster, key=lambda x: counts[x])
        else:
            # Use most common UMI as representative
            key = max(cluster, key=lambda x: counts[x])
            
        # Map all UMIs in cluster to representative
        for item in cluster:
            mapping[item] = key
            collapsed_counts[key] = collapsed_counts.get(key, 0) + counts[item]
    
    return collapsed_counts, mapping


class UMI:
    """Class for handling UMI clustering and correction"""
    
    def __init__(self, mismatches: int = 2, method: str = "adjacency", correct: bool = True):
        """Initialize UMI clusterer
        
        Args:
            mismatches: Number of mismatches allowed for clustering
            method: Clustering method to use (cluster, adjacency, or directional)
            correct: Whether to perform clustering (default: True)
        """
        self.clusterer = UMIClusterer(cluster_method=method)
        self.mismatches = mismatches
        self.method = method  # Store the method name
        self._counts: Dict[bytes, int] = {}
        self._mapping: Dict[bytes, bytes] = {}
        self._corrected_counts: Optional[Dict[bytes, int]] = None
        self.correct = correct
    
    def consume(self, umi: Union[str, bytes], n: int = 1) -> None:
        """Add a UMI to the counts dictionary
        
        Args:
            umi: UMI sequence to add
            n: Number of times to add the UMI (default: 1)
        """
        if isinstance(umi, str):
            umi = umi.encode('ascii')
        self._counts[umi] = self._counts.get(umi, 0) + n
        self._corrected_counts = None  # Reset corrected counts
    
    def create_mapping(self, allowed_list: Optional[List[Union[str, bytes]]] = None) -> None:
        """Create mapping between original and corrected barcodes
        
        Args:
            allowed_list: Optional list of allowed UMIs. If provided, highest priority UMI in cluster will be used as representative
        """
        if not self._counts:
            return
            
        # Convert allowed_list to bytes if provided
        if allowed_list:
            allowed_list = [x.encode('ascii') if isinstance(x, str) else x for x in allowed_list]
            
        # If correction disabled, each barcode maps to itself
        if not self.correct:
            self._mapping = {bc: bc for bc in self._counts.keys()}
            self._corrected_counts = self._counts.copy()
            return
            
        # Use private clustering function
        self._corrected_counts, self._mapping = _cluster_umis(
            self._counts, 
            self.mismatches, 
            self.method,  # Use stored method name
            allowed_list
        )
    
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
            # If no mapping exists, use original counts
            if not self._mapping:
                self._corrected_counts = self._counts.copy()
            else:
                # Create corrected counts from mapping
                self._corrected_counts = {}
                for orig_bc, count in self._counts.items():
                    corrected = self._mapping.get(orig_bc, orig_bc)
                    self._corrected_counts[corrected] = self._corrected_counts.get(corrected, 0) + count
        return self._corrected_counts
    
    @classmethod
    def from_csv(cls, filepath: Union[str, Path], column: str, 
                 mismatches: int = 2, method: str = "adjacency",
                 sep: str = ",", correct: bool = True,
                 count_column: Optional[str] = None,
                 scale: Optional[float] = None) -> 'UMI':
        """Create UMI object from CSV file
        
        Args:
            filepath: Path to CSV file
            column: Column containing UMIs
            mismatches: Number of mismatches allowed for clustering
            method: Clustering method to use
            sep: CSV separator
            correct: If True, perform clustering. If False, assume data is pre-clustered.
            count_column: Optional column containing pre-counted values
            scale: Optional scale factor for normalized values
            
        Returns:
            UMI object
        """
        umi = cls(mismatches=mismatches, method=method, correct=correct)
        
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f, delimiter=sep)
            if column not in reader.fieldnames:
                raise ValueError(f"Column {column} not found in CSV file")
            if count_column and count_column not in reader.fieldnames:
                raise ValueError(f"Count column {count_column} not found in CSV file")
            
            for row in reader:
                value = row[column]
                if not value:  # Skip empty values
                    continue
                    
                if count_column:
                    try:
                        count = float(row[count_column])
                        if scale:
                            count = int(round(count * scale))
                        else:
                            count = int(round(count))
                    except (ValueError, TypeError):
                        print(f"Warning: Invalid count value '{row[count_column]}' for {value}", file=sys.stderr)
                        continue
                    
                    # Add the value count times
                    umi.consume(value, count)
                else:
                    umi.consume(value)
        
        if correct:
            umi.create_mapping()
        else:
            # For pre-clustered data, use each UMI as its own cluster
            umi._mapping = {bc: bc for bc in umi._counts.keys()}
            umi._corrected_counts = umi._counts.copy()
        
        return umi

class UmiCollection:
    """Class for managing collections of UMI objects across samples"""
    
    def __init__(self, umis: Dict[str, UMI] = None):
        """Initialize UmiCollection
        
        Args:
            umis: Dictionary mapping sample names to UMI objects
        """
        self.umis = umis or {}
        self.logger = logging.getLogger(__name__)
    
    def collapse_umis(self, mismatches: int = 2, method: str = "adjacency", 
                     allowed_list: Optional[List[Union[str, bytes]]] = None) -> 'UmiCollection':
        """Collapse UMIs across all samples using clustering
        
        Args:
            mismatches: Number of mismatches allowed for clustering
            method: Clustering method to use
            allowed_list: Optional list of allowed UMIs. If provided, highest priority UMI in cluster will be used as representative
            
        Returns:
            New UmiCollection with collapsed UMIs
        """
        # Combine all UMIs from all samples
        combined_counts = {}
        for umi_obj in self.umis.values():
            for bc, count in umi_obj._counts.items():
                combined_counts[bc] = combined_counts.get(bc, 0) + count
        
        # Convert allowed_list to bytes if provided
        if allowed_list:
            allowed_list = [x.encode('ascii') if isinstance(x, str) else x for x in allowed_list]
        
        # Use private clustering function
        collapsed_counts, mapping = _cluster_umis(
            combined_counts,
            mismatches,
            method,
            allowed_list
        )
        
        # Create new collection with collapsed counts
        new_umis = {}
        for sample_name, umi_obj in self.umis.items():
            new_umi = UMI(mismatches=mismatches, method=method, correct=False)
            new_umi._mapping = mapping
            new_umi._counts = umi_obj._counts.copy()
            
            # Update corrected counts using the mapping
            new_umi._corrected_counts = {}
            for orig_bc, count in umi_obj._counts.items():
                corrected = mapping.get(orig_bc, orig_bc)
                new_umi._corrected_counts[corrected] = new_umi._corrected_counts.get(corrected, 0) + count
                
            new_umis[sample_name] = new_umi
            
        return UmiCollection(new_umis)
    
    @classmethod
    def from_csvs(cls, filepaths: List[Union[str, Path]], 
                  column: str,
                  sample_names: Optional[List[str]] = None,
                  mismatches: int = 0,
                  method: str = "adjacency",
                  sep: str = ",",
                  count_column: Optional[str] = None) -> 'UmiCollection':
        """Create UmiCollection from multiple CSV files
        
        Args:
            filepaths: List of paths to CSV files
            column: Column containing UMIs
            sample_names: Optional list of sample names. If not provided, uses basenames
            mismatches: Number of mismatches allowed for clustering
            method: Clustering method to use
            sep: CSV separator
            count_column: Optional column containing counts
            
        Returns:
            UmiCollection object
        """
        logger = logging.getLogger(__name__)
        
        if sample_names and len(sample_names) != len(filepaths):
            raise ValueError("Number of sample names must match number of files")
            
        umis = {}
        for i, filepath in enumerate(filepaths):
            # Use provided sample name or basename
            sample_name = sample_names[i] if sample_names else Path(filepath).stem
            logger.info(f"Processing file {i+1}/{len(filepaths)}: {filepath} as sample {sample_name}")
            
            # Create UMI object for this file
            umi = UMI(mismatches=mismatches, method=method, correct=False)
            
            # Read CSV and add counts
            with open(filepath, 'r') as f:
                reader = csv.DictReader(f, delimiter=sep)
                if column not in reader.fieldnames:
                    raise ValueError(f"Column {column} not found in {filepath}")
                if count_column and count_column not in reader.fieldnames:
                    raise ValueError(f"Count column {count_column} not found in {filepath}")
                
                for row in reader:
                    value = row[column]
                    if not value:  # Skip empty values
                        continue
                        
                    if count_column:
                        try:
                            count = int(float(row[count_column]))
                            umi.consume(value, count)
                        except (ValueError, TypeError):
                            logger.warning(f"Invalid count value '{row[count_column]}' for {value}")
                            continue
                    else:
                        umi.consume(value)
            
            umis[sample_name] = umi
            
        return cls(umis)
    
    @classmethod
    def from_df(cls, df: 'pd.DataFrame',
                sample_col: str,
                umi_col: str,
                count_col: Optional[str] = None) -> 'UmiCollection':
        """Create UmiCollection from pandas DataFrame
        
        Args:
            df: DataFrame containing UMI data
            sample_col: Column containing sample names
            umi_col: Column containing UMIs
            count_col: Optional column containing counts
            
        Returns:
            UmiCollection object
        """
        logger = logging.getLogger(__name__)
        umis = {}
        
        for sample in df[sample_col].unique():
            logger.info(f"Processing sample: {sample}")
            sample_df = df[df[sample_col] == sample]
            umi = UMI(mismatches=0, correct=False)
            
            if count_col:
                for _, row in sample_df.iterrows():
                    umi.consume(row[umi_col], int(row[count_col]))
            else:
                for umi_seq in sample_df[umi_col]:
                    umi.consume(umi_seq)
                    
            umis[sample] = umi
            
        return cls(umis)
    
    def to_df(self, format: str = 'wide') -> 'pd.DataFrame':
        """Convert UmiCollection to pandas DataFrame
        
        Args:
            format: Output format, either 'wide' or 'long'
                  - wide: Each sample is a column, rows are UMIs
                  - long: Three columns: sample, umi, count
        
        Returns:
            DataFrame in specified format
        """
        import pandas as pd
        
        if format not in ['wide', 'long']:
            raise ValueError("Format must be either 'wide' or 'long'")
            
    
        rows = []
        for sample, umi in self.umis.items():
            for bc, count in umi.corrected_counts.items():
                rows.append({
                    'sample': sample,
                    'umi': bc.decode('ascii'),
                    'count': count
                })
        long = pd.DataFrame(rows)
        if format == 'wide':
            return pd.pivot_table(long,
                                    values='count',
                                    index='umi',
                                    columns='sample',
                                    fill_value=0)
                
        return long
    
    def write(self, path: Union[str, Path], sep: str = ",", format: str = 'wide') -> None:
        """Write UmiCollection to CSV file
        
        Args:
            path: Path to output CSV file
            sep: CSV separator
            format: Output format, either 'wide' or 'long'
        """
        self.logger.info(f"Writing data to {path} in {format} format")
        df = self.to_df(format=format)
        df.to_csv(path, sep=sep, index=(format == 'wide'))
        self.logger.info(f"Successfully wrote {len(df)} rows to {path}") 