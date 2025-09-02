"""Tests for UMI clustering and correction"""

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
from outerspace.umi import UMI, UmiCollection
from pathlib import Path
import tempfile
import csv
import os
import pandas as pd
import numpy as np


def test_basic_umi_merging():
    """Test basic UMI merging with single mismatch"""
    umi = UMI(mismatches=2)

    # Add UMIs that should be merged
    umi.consume("ATCG")
    umi.consume("ATCC")  # One mismatch from ATCG
    umi.consume("GCTA")  # Different UMI

    umi.create_mapping()

    # Check that similar UMIs are merged
    assert umi["ATCG"] == umi["ATCC"]  # Should be merged
    assert umi["GCTA"] != umi["ATCG"]  # Should not be merged

    # Check counts
    counts = umi.corrected_counts
    assert len(counts) == 2  # Should have 2 unique clusters
    assert sum(counts.values()) == 3  # Total count should be 3


def test_pre_clustered_data():
    """Test loading pre-clustered data"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        writer = csv.writer(f)
        writer.writerow(["barcode"])
        writer.writerow(["ATCG"])
        writer.writerow(["ATCC"])
        writer.writerow(["GCTA"])
        temp_file = f.name

    try:
        # Load with correction disabled
        umi = UMI.from_csv(temp_file, column="barcode", correct=False)

        # Check that UMIs are not merged
        assert umi["ATCG"] != umi["ATCC"]  # Should not be merged
        assert umi["GCTA"] != umi["ATCG"]  # Should not be merged

        # Check counts - should be same as original
        counts = umi.corrected_counts
        assert len(counts) == 3  # Should have 3 unique clusters
        assert sum(counts.values()) == 3  # Total count should be 3

        # Verify mapping is identity
        assert all(umi[bc] == bc for bc in counts.keys())

        # Compare with corrected version
        umi_corrected = UMI.from_csv(temp_file, column="barcode", correct=True)
        assert len(umi_corrected.corrected_counts) < len(
            counts
        )  # Should have fewer clusters

    finally:
        os.unlink(temp_file)


def test_csv_loading():
    """Test loading UMIs from CSV file"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        writer = csv.writer(f)
        writer.writerow(["barcode"])
        writer.writerow(["ATCG"])
        writer.writerow(["ATCC"])
        writer.writerow(["GCTA"])
        temp_file = f.name

    try:
        # Test single file loading
        umi = UMI.from_csv(temp_file, column="barcode", mismatches=0, correct=False)
        counts = umi.corrected_counts
        assert len(counts) == 3  # Three unique barcodes
        assert sum(counts.values()) == 3  # One count each

        # Test multiple file loading
        temp_dir = tempfile.mkdtemp()
        try:
            # Create multiple CSV files
            for i in range(2):
                with open(os.path.join(temp_dir, f"test_{i}.csv"), "w") as f:
                    writer = csv.writer(f)
                    writer.writerow(["barcode"])
                    writer.writerow(["ATCG"])
                    writer.writerow(["ATCC"])

            # Use UmiCollection for multiple files
            collection = UmiCollection.from_csvs(
                [os.path.join(temp_dir, f"test_{i}.csv") for i in range(2)],
                column="barcode",
                sample_names=["sample1", "sample2"],
            )
            assert len(collection.umis) == 2
            for sample in collection.umis.values():
                counts = sample.corrected_counts
                assert len(counts) == 2  # Two unique barcodes per sample
                assert sum(counts.values()) == 2  # One count each

        finally:
            # Cleanup
            for file in Path(temp_dir).glob("*.csv"):
                file.unlink()
            os.rmdir(temp_dir)

    finally:
        os.unlink(temp_file)


def test_empty_input():
    """Test handling of empty input"""
    umi = UMI()
    umi.create_mapping()
    assert len(umi.corrected_counts) == 0

    # Test empty CSV
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        writer = csv.writer(f)
        writer.writerow(["barcode"])
        temp_file = f.name

    try:
        umi = UMI.from_csv(temp_file, column="barcode")
        assert len(umi.corrected_counts) == 0
    finally:
        os.unlink(temp_file)


def test_invalid_input():
    """Test handling of invalid input"""
    umi = UMI()

    # Test invalid characters
    with pytest.raises(UnicodeEncodeError):
        umi.consume("ATCG\u2603")  # Snowman character

    # Test invalid CSV column
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        writer = csv.writer(f)
        writer.writerow(["wrong_column"])
        writer.writerow(["ATCG"])
        temp_file = f.name

    try:
        with pytest.raises(ValueError):
            UMI.from_csv(temp_file, column="barcode")
    finally:
        os.unlink(temp_file)


def test_umi_collection_from_csvs():
    """Test creating UmiCollection from multiple CSV files"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test CSV files with different sample data
        files = []
        for i in range(2):
            file_path = os.path.join(temp_dir, f"sample_{i}.csv")
            with open(file_path, "w") as f:
                writer = csv.writer(f)
                writer.writerow(["barcode", "count"])
                writer.writerow(["ATCG", "2"])
                writer.writerow(["ATCC", "1"])
                writer.writerow(["GCTA", "3"])
            files.append(file_path)

        # Create collection from files
        collection = UmiCollection.from_csvs(
            files,
            column="barcode",
            count_column="count",
            sample_names=["sample1", "sample2"],
        )

        # Check that we have the right number of samples
        assert len(collection.umis) == 2

        # Check counts for each sample
        for sample in ["sample1", "sample2"]:
            counts = collection.umis[sample].corrected_counts
            assert len(counts) == 3  # Three unique barcodes
            assert sum(counts.values()) == 6  # Total count of 6 (2+1+3)


def test_umi_collection_from_df():
    """Test creating UmiCollection from pandas DataFrame"""
    # Create test DataFrame
    df = pd.DataFrame(
        {
            "sample": ["A", "A", "A", "B", "B", "B"],
            "barcode": ["ATCG", "ATCC", "GCTA", "ATCG", "ATCC", "GCTA"],
            "count": [2, 1, 3, 1, 2, 2],
        }
    )

    # Create collection from DataFrame
    collection = UmiCollection.from_df(
        df, sample_col="sample", umi_col="barcode", count_col="count"
    )

    # Check that we have the right number of samples
    assert len(collection.umis) == 2

    # Check counts for each sample
    for sample in ["A", "B"]:
        counts = collection.umis[sample].corrected_counts
        assert len(counts) == 3  # Three unique barcodes
        if sample == "A":
            assert sum(counts.values()) == 6  # 2+1+3
        else:
            assert sum(counts.values()) == 5  # 1+2+2


def test_umi_collection_to_df_wide():
    """Test converting UmiCollection to wide format DataFrame"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test CSV files
        files = []
        for i in range(2):
            file_path = os.path.join(temp_dir, f"sample_{i}.csv")
            with open(file_path, "w") as f:
                writer = csv.writer(f)
                writer.writerow(["barcode", "count"])
                writer.writerow(["ATCG", "2"])
                writer.writerow(["ATCC", "1"])
                writer.writerow(["GCTA", "3"])
            files.append(file_path)

        # Create collection
        collection = UmiCollection.from_csvs(
            files,
            column="barcode",
            count_column="count",
            sample_names=["sample1", "sample2"],
        )

        # Convert to wide format DataFrame
        df = collection.to_df(format="wide")

        # Check DataFrame structure
        assert isinstance(df, pd.DataFrame)
        assert df.shape == (3, 2)  # 3 barcodes, 2 samples
        assert list(df.columns) == ["sample1", "sample2"]
        assert set(df.index) == {
            "ATCG",
            "ATCC",
            "GCTA",
        }  # Check set of indices instead of order

        # Check values
        assert df.loc["ATCG", "sample1"] == 2
        assert df.loc["ATCC", "sample1"] == 1
        assert df.loc["GCTA", "sample1"] == 3


def test_umi_collection_to_df_long():
    """Test converting UmiCollection to long format DataFrame"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test CSV files
        files = []
        for i in range(2):
            file_path = os.path.join(temp_dir, f"sample_{i}.csv")
            with open(file_path, "w") as f:
                writer = csv.writer(f)
                writer.writerow(["barcode", "count"])
                writer.writerow(["ATCG", "2"])
                writer.writerow(["ATCC", "1"])
                writer.writerow(["GCTA", "3"])
            files.append(file_path)

        # Create collection
        collection = UmiCollection.from_csvs(
            files,
            column="barcode",
            count_column="count",
            sample_names=["sample1", "sample2"],
        )

        # Convert to long format DataFrame
        df = collection.to_df(format="long")

        # Check DataFrame structure
        assert isinstance(df, pd.DataFrame)
        assert df.shape == (6, 3)  # 6 rows (3 barcodes × 2 samples), 3 columns
        assert list(df.columns) == ["sample", "umi", "count"]

        # Check values
        sample1_data = df[df["sample"] == "sample1"]
        assert len(sample1_data) == 3
        assert sample1_data[sample1_data["umi"] == "ATCG"]["count"].iloc[0] == 2
        assert sample1_data[sample1_data["umi"] == "ATCC"]["count"].iloc[0] == 1
        assert sample1_data[sample1_data["umi"] == "GCTA"]["count"].iloc[0] == 3


def test_umi_collection_write():
    """Test writing UmiCollection to CSV file"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test CSV files
        files = []
        for i in range(2):
            file_path = os.path.join(temp_dir, f"sample_{i}.csv")
            with open(file_path, "w") as f:
                writer = csv.writer(f)
                writer.writerow(["barcode", "count"])
                writer.writerow(["ATCG", "2"])
                writer.writerow(["ATCC", "1"])
                writer.writerow(["GCTA", "3"])
            files.append(file_path)

        # Create collection
        collection = UmiCollection.from_csvs(
            files,
            column="barcode",
            count_column="count",
            sample_names=["sample1", "sample2"],
        )

        # Test writing in wide format
        wide_output = os.path.join(temp_dir, "wide_output.csv")
        collection.write(wide_output, format="wide")
        assert os.path.exists(wide_output)

        # Read back and verify
        wide_df = pd.read_csv(wide_output, index_col=0)
        assert wide_df.shape == (3, 2)
        assert list(wide_df.columns) == ["sample1", "sample2"]

        # Test writing in long format
        long_output = os.path.join(temp_dir, "long_output.csv")
        collection.write(long_output, format="long")
        assert os.path.exists(long_output)

        # Read back and verify
        long_df = pd.read_csv(long_output)
        assert long_df.shape == (6, 3)
        assert list(long_df.columns) == ["sample", "umi", "count"]


def test_umi_collection_invalid_format():
    """Test that UmiCollection handles invalid format specifications"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test CSV file
        file_path = os.path.join(temp_dir, "sample.csv")
        with open(file_path, "w") as f:
            writer = csv.writer(f)
            writer.writerow(["barcode", "count"])
            writer.writerow(["ATCG", "2"])

        # Create collection
        collection = UmiCollection.from_csvs(
            [file_path],
            column="barcode",
            count_column="count",
            sample_names=["sample1"],
        )

        # Test invalid format in to_df
        with pytest.raises(ValueError, match="Format must be either 'wide' or 'long'"):
            collection.to_df(format="invalid")

        # Test invalid format in write
        with pytest.raises(ValueError, match="Format must be either 'wide' or 'long'"):
            collection.write(os.path.join(temp_dir, "output.csv"), format="invalid")


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.