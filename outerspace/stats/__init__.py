"""Statistics module for UMI analysis"""

from .base import BaseStatistic, BasePairwiseStatistic, BaseDifferentialStatistic

# Import statistics as they are implemented
from .single import (
    GiniCoefficient,
    ShannonDiversity,
    SimpsonDiversity,
    UMIRecoveryRate,
    UMIEfficiencyRate,
    UMIErrorRate,
    UMIRedundancy,
)

from .multi import (
    JaccardSimilarity,
    BrayCurtisDissimilarity,
    SpearmanCorrelation,
    FoldChange,
)

from .differential import (
    SingleSampleDifferentialAbundance,
    MannWhitneyDifferentialAbundance,
    PairedTTestDifferentialAbundance,
)


__all__ = [
    "BaseStatistic",
    "BasePairwiseStatistic",
    "BaseDifferentialStatistic",
    "GiniCoefficient",
    "ShannonDiversity",
    "SimpsonDiversity",
    "UMIRecoveryRate",
    "UMIEfficiencyRate",
    "UMIErrorRate",
    "UMIRedundancy",
    "JaccardSimilarity",
    "BrayCurtisDissimilarity",
    "SpearmanCorrelation",
    "FoldChange",
    "SingleSampleDifferentialAbundance",
    "MannWhitneyDifferentialAbundance",
    "PairedTTestDifferentialAbundance"
    # Add other classes as they are implemented
]
