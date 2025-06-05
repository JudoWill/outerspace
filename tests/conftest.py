"""Shared pytest fixtures for OUTERSPACE tests"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import os
import pytest
import tempfile
import shutil

@pytest.fixture
def temp_workspace():
    """Create a temporary workspace with test data structure"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create directory structure
        os.makedirs(os.path.join(temp_dir, 'reads/'))
        os.makedirs(os.path.join(temp_dir, 'findseq'))
        os.makedirs(os.path.join(temp_dir, 'collapse'))
        os.makedirs(os.path.join(temp_dir, 'count'))
        os.makedirs(os.path.join(temp_dir, 'gini'))
        
        # Create symbolic link to workflow directory
        workflow_link = os.path.join(temp_dir, 'workflow')
        os.symlink(os.path.abspath('workflow'), workflow_link)
        
        # Copy the real data files to temp directory
        data_files = [
            ('tests/data/409-4_S1_L002_R1_001.fastq.gz', 'reads/409-4_S1_L002_R1_001.fastq.gz'),
            ('tests/data/409-4_S1_L002_R2_001.fastq.gz', 'reads/409-4_S1_L002_R2_001.fastq.gz'),
            ('tests/data/2-G1L9-M1_S9_L001_R1_001.fastq.gz', 'reads/2-G1L9-M1_S9_L001_R1_001.fastq.gz'),
            ('tests/data/2-G1L9-M1_S9_L001_R2_001.fastq.gz', 'reads/2-G1L9-M1_S9_L001_R2_001.fastq.gz'),
            ('tests/data/2-G1L9-M2_S12_L001_R1_001.fastq.gz', 'reads/2-G1L9-M2_S12_L001_R1_001.fastq.gz'),
            ('tests/data/2-G1L9-M2_S12_L001_R2_001.fastq.gz', 'reads/2-G1L9-M2_S12_L001_R2_001.fastq.gz'),
            ('tests/data/library_protospacers.txt', 'library_protospacers.txt')
        ]
        
        for src, dst in data_files:
            src_path = src
            dst_path = os.path.join(temp_dir, dst)
            if os.path.exists(src_path):
                shutil.copy2(src_path, dst_path)
            else:
                # Create empty files for testing if real files don't exist
                os.makedirs(os.path.dirname(dst_path), exist_ok=True)
                open(dst_path, 'w').close()
        
        # Copy the config files
        config_files = [
            ('tests/test_cli/pipeline_data/samplesheet_config.yaml', 'samplesheet_config.yaml'),
            ('tests/test_cli/pipeline_data/samplesheet.csv', 'samplesheet.csv'),
            ('tests/test_cli/pipeline_data/directory_config.yaml', 'directory_config.yaml'),
            ('tests/test_cli/pipeline_data/grnaquery.toml', 'grnaquery.toml')
        ]
        
        for src, dst in config_files:
            src_path = src
            dst_path = os.path.join(temp_dir, dst)
            if os.path.exists(src_path):
                shutil.copy2(src_path, dst_path)
            else:
                raise FileNotFoundError(f"Could not find config file: {src_path} for testing")

        yield temp_dir 