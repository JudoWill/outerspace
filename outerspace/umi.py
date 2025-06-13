"""Module for handling UMI clustering and correction"""

from typing import Dict, List, Optional, Union
from collections import Counter
from umi_tools import UMIClusterer
import csv
import os
from pathlib import Path
import sys
import logging


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
    
    def create_mapping(self) -> None:
        """Create mapping between original and corrected barcodes"""
        if not self._counts:
            return
            
        # If no mismatches allowed or correction disabled, each barcode maps to itself
        if self.mismatches == 0 or not self.correct:
            self._mapping = {bc: bc for bc in self._counts.keys()}
            self._corrected_counts = self._counts.copy()
            return
            
        # Create clusters
        clusters = self.clusterer(self._counts, self.mismatches)
        
        # Create mapping
        self._mapping = {}
        for cluster in clusters:
            key = cluster[0]  # Use first barcode in cluster as representative
            for item in cluster:
                self._mapping[item] = key
        
        # Update corrected counts
        self._corrected_counts = {}
        for orig_bc, count in self._counts.items():
            corrected = self._mapping.get(orig_bc, orig_bc)
            self._corrected_counts[corrected] = self._corrected_counts.get(corrected, 0) + count
    
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