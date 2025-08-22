"""Differential abundance calculations for UMI analysis.

This module provides classes for calculating differential abundance statistics
between groups of UMI libraries, including various statistical tests and
effect size measures.
"""

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats

from outerspace.stats.base import BaseDifferentialStatistic
from outerspace.stats.utils import split_counts_by_allowed_list
from outerspace.umi import UMI

# Set up logging
logger = logging.getLogger(__name__)

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


class BaseDifferentialAbundance(BaseDifferentialStatistic):
    """Base class for differential abundance calculations.

    This class provides common functionality for differential abundance analysis
    including data preparation and basic statistical calculations.
    """

    def __init__(
        self,
        umis: List[UMI],
        groups: List[str],
        use_corrected: bool = True,
        allowed_list: Optional[List[str]] = None,
    ) -> None:
        """Initialize differential abundance calculator.

        Parameters
        ----------
        umis : List[UMI]
            List of UMI objects
        groups : List[str]
            List of group labels corresponding to UMI objects
        use_corrected : bool, default=True
            If True, use corrected counts. If False, use original counts.
        allowed_list : Optional[List[str]], default=None
            Optional list of allowed UMIs. If provided, only these UMIs will be considered.
        """
        super().__init__(umis, groups)
        self.use_corrected = use_corrected
        self.allowed_list = allowed_list

    @property
    def _allowed_list(self) -> Optional[List[bytes]]:
        """Get allowed list in bytes format.

        Returns
        -------
        Optional[List[bytes]]
            Allowed list converted to bytes, or None if not provided
        """
        if self.allowed_list:
            return [umi.encode("ascii") for umi in self.allowed_list]
        return None

    @staticmethod
    def _calculate_log2_fold_change(mean1: float, mean2: float) -> float:
        """Calculate log2 fold change between two means.

        Parameters
        ----------
        mean1 : float
            First mean value
        mean2 : float
            Second mean value

        Returns
        -------
        float
            Log2 fold change

        Notes
        -----
        If mean1 is 0 and mean2 > 0, returns infinity.
        If both means are 0, returns 0.
        """
        if mean1 == 0:
            return float("inf") if mean2 > 0 else 0.0
        return np.log2(mean2 / mean1)

    def _prepare_dataframes(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Prepare DataFrames for differential abundance calculation.

        Returns
        -------
        Tuple[pd.DataFrame, pd.DataFrame]
            Tuple of (group1_df, group2_df) containing counts for each group

        Raises
        ------
        ValueError
            If there are not exactly 2 groups
        """
        # Create DataFrame of counts
        all_umis = set()
        for umi in self.umis:
            counts = umi.corrected_counts if self.use_corrected else umi._counts
            if self.allowed_list:
                counts, _, _ = split_counts_by_allowed_list(
                    counts, self._allowed_list, add_missing=True
                )
            all_umis.update(counts.keys())

        # Create DataFrames for each group
        group1_df = pd.DataFrame(index=list(all_umis))
        group2_df = pd.DataFrame(index=list(all_umis))

        # Get unique groups
        unique_groups = sorted(set(self.groups))
        if len(unique_groups) != 2:
            raise ValueError(
                "Differential abundance calculation requires exactly 2 groups"
            )
        group1, group2 = unique_groups

        # Add samples to appropriate group DataFrame
        for i, (umi, group) in enumerate(zip(self.umis, self.groups)):
            counts = umi.corrected_counts if self.use_corrected else umi._counts
            if self.allowed_list:
                counts, _, _ = split_counts_by_allowed_list(
                    counts, self._allowed_list, add_missing=True
                )

            if group == group1:
                group1_df[f"sample_{i}"] = pd.Series(counts)
            else:
                group2_df[f"sample_{i}"] = pd.Series(counts)

        logger.debug(
            f"Prepared DataFrames: group1 ({group1}) with {group1_df.shape[1]} samples, group2 ({group2}) with {group2_df.shape[1]} samples"
        )
        return group1_df, group2_df


class SingleSampleDifferentialAbundance(BaseDifferentialAbundance):
    """Calculate differential abundance scores between single samples.

    This class provides differential abundance analysis for single-sample
    comparisons, useful for pilot studies or when replicates are not available.
    """

    @staticmethod
    def calculate_differential_abundance_single(
        group1_df: pd.DataFrame, group2_df: pd.DataFrame
    ) -> pd.DataFrame:
        """Calculate differential abundance scores between single samples.

        Parameters
        ----------
        group1_df : pd.DataFrame
            DataFrame containing counts for first group
        group2_df : pd.DataFrame
            DataFrame containing counts for second group

        Returns
        -------
        pd.DataFrame
            DataFrame containing differential abundance metrics for each UMI

        Raises
        ------
        ValueError
            If DataFrames are empty or don't have exactly one sample each
        """
        # Validate DataFrames
        if group1_df.empty or group2_df.empty:
            raise ValueError("Empty DataFrame provided")
        if group1_df.shape[1] != 1 or group2_df.shape[1] != 1:
            raise ValueError(
                "Single-sample comparison requires exactly one sample per group"
            )

        # Calculate metrics for each UMI
        results = []
        for umi in group1_df.index:
            # Get counts for each group
            group1_count = group1_df.loc[umi].values[0]
            group2_count = group2_df.loc[umi].values[0]

            # Calculate log2 fold change
            lfc = BaseDifferentialAbundance._calculate_log2_fold_change(
                group1_count, group2_count
            )

            # Calculate effect size (absolute difference since we can't calculate std with one sample)
            effect_size = abs(group2_count - group1_count)

            results.append(
                {
                    "umi": umi,
                    "log2_fold_change": lfc,
                    "effect_size": effect_size,
                    "mean_group1": group1_count,
                    "mean_group2": group2_count,
                    "std_group1": 0.0,  # No std with one sample
                    "std_group2": 0.0,  # No std with one sample
                }
            )

        return pd.DataFrame(results).set_index("umi")

    def run(self) -> pd.DataFrame:
        """Calculate differential abundance scores between single samples.

        Returns
        -------
        pd.DataFrame
            DataFrame containing differential abundance metrics for each UMI
        """
        group1_df, group2_df = self._prepare_dataframes()
        result = self.calculate_differential_abundance_single(group1_df, group2_df)
        logger.info(
            f"Calculated differential abundance for {len(result)} UMIs between single samples"
        )
        return result


class MannWhitneyDifferentialAbundance(BaseDifferentialAbundance):
    """Calculate differential abundance scores using Mann-Whitney U test.

    This class provides non-parametric differential abundance analysis using
    the Mann-Whitney U test, which is robust to non-normal distributions.
    """

    @staticmethod
    def calculate_differential_abundance(
        group1_df: pd.DataFrame, group2_df: pd.DataFrame
    ) -> pd.DataFrame:
        """Calculate differential abundance scores using Mann-Whitney U test.

        Parameters
        ----------
        group1_df : pd.DataFrame
            DataFrame containing counts for first group
        group2_df : pd.DataFrame
            DataFrame containing counts for second group

        Returns
        -------
        pd.DataFrame
            DataFrame containing differential abundance metrics for each UMI

        Raises
        ------
        ValueError
            If DataFrames are empty
        """
        # Validate DataFrames
        if group1_df.empty or group2_df.empty:
            raise ValueError("Empty DataFrame provided")

        # Calculate metrics for each UMI
        results = []
        for umi in group1_df.index:
            # Get counts for each group
            group1_counts = group1_df.loc[umi].values
            group2_counts = group2_df.loc[umi].values

            # Calculate means and standard deviations
            mean1 = np.mean(group1_counts)
            mean2 = np.mean(group2_counts)
            std1 = np.std(group1_counts)
            std2 = np.std(group2_counts)

            # Calculate log2 fold change
            lfc = BaseDifferentialAbundance._calculate_log2_fold_change(mean1, mean2)

            # Calculate effect size (Cohen's d)
            denominator = np.sqrt((std1**2 + std2**2) / 2)
            if denominator == 0:
                effect_size = 0.0
            else:
                effect_size = (mean2 - mean1) / denominator

            # Calculate p-value using Mann-Whitney U test
            _, p_value = stats.mannwhitneyu(
                group1_counts, group2_counts, alternative="two-sided"
            )

            results.append(
                {
                    "umi": umi,
                    "log2_fold_change": lfc,
                    "effect_size": effect_size,
                    "p_value": p_value,
                    "mean_group1": mean1,
                    "mean_group2": mean2,
                    "std_group1": std1,
                    "std_group2": std2,
                }
            )

        return pd.DataFrame(results).set_index("umi")

    def run(self) -> pd.DataFrame:
        """Calculate differential abundance scores using Mann-Whitney U test.

        Returns
        -------
        pd.DataFrame
            DataFrame containing differential abundance metrics for each UMI
        """
        group1_df, group2_df = self._prepare_dataframes()
        result = self.calculate_differential_abundance(group1_df, group2_df)
        logger.info(
            f"Calculated Mann-Whitney differential abundance for {len(result)} UMIs"
        )
        return result


class PairedTTestDifferentialAbundance(BaseDifferentialAbundance):
    """Calculate differential abundance scores using paired t-test.

    This class provides parametric differential abundance analysis using
    paired t-tests, suitable for before/after or matched-pair experimental designs.
    """

    @staticmethod
    def calculate_differential_abundance_paired(
        pre_df: pd.DataFrame, post_df: pd.DataFrame
    ) -> pd.DataFrame:
        """Calculate differential abundance scores using paired t-test.

        Parameters
        ----------
        pre_df : pd.DataFrame
            DataFrame containing pre-treatment counts
        post_df : pd.DataFrame
            DataFrame containing post-treatment counts

        Returns
        -------
        pd.DataFrame
            DataFrame containing differential abundance metrics for each UMI

        Raises
        ------
        ValueError
            If DataFrames are empty or have different numbers of samples
        """
        # Validate DataFrames
        if pre_df.empty or post_df.empty:
            raise ValueError("Empty DataFrame provided")
        if pre_df.shape[1] != post_df.shape[1]:
            raise ValueError("Pre and post DataFrames must have same number of samples")

        # Calculate metrics for each UMI
        results = []
        for umi in pre_df.index:
            # Get counts for each group
            pre_counts = pre_df.loc[umi].values
            post_counts = post_df.loc[umi].values

            # Calculate means and standard deviations
            mean_pre = np.mean(pre_counts)
            mean_post = np.mean(post_counts)
            std_pre = np.std(pre_counts)
            std_post = np.std(post_counts)

            # Calculate log2 fold change
            lfc = BaseDifferentialAbundance._calculate_log2_fold_change(
                mean_pre, mean_post
            )

            # Calculate effect size (Cohen's d for paired samples)
            differences = post_counts - pre_counts
            if np.std(differences) > 0:
                effect_size = np.mean(differences) / np.std(differences)
            else:
                effect_size = np.mean(differences)

            # Calculate p-value using paired t-test
            # If all values are identical, p-value should be 1.0
            if np.all(pre_counts == post_counts):
                p_value = 1.0
            else:
                _, p_value = stats.ttest_rel(post_counts, pre_counts)

            results.append(
                {
                    "umi": umi,
                    "log2_fold_change": lfc,
                    "effect_size": effect_size,
                    "p_value": p_value,
                    "mean_pre": mean_pre,
                    "mean_post": mean_post,
                    "std_pre": std_pre,
                    "std_post": std_post,
                }
            )

        return pd.DataFrame(results).set_index("umi")

    def run(self) -> pd.DataFrame:
        """Calculate differential abundance scores using paired t-test.

        Returns
        -------
        pd.DataFrame
            DataFrame containing differential abundance metrics for each UMI
        """
        group1_df, group2_df = self._prepare_dataframes()
        result = self.calculate_differential_abundance_paired(group1_df, group2_df)
        logger.info(
            f"Calculated paired t-test differential abundance for {len(result)} UMIs"
        )
        return result


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
