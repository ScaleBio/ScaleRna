import numpy as np
from dataclasses import dataclass


@dataclass(frozen=True)
class OutlierOptions:
    """
    Options for flagging / filtering cell outliers based on QC metrics
    """

    filter_outliers: bool
    reads_mads: float
    passing_mads: float
    mito_mads: float
    mito_min_thres: float = 0.05


def coverage_equation(x, c, n):
    """
    Equation for expected coverage at a given read-number.
    """
    return c / x - 1 + np.exp(-n / x)


def estLibSize(reads: int, umis: int) -> float:
    """
    Function to estimate library size(number of unique molecules)
    at infinite sequencing depth
    Based on the method from picard tools
    Args:
        reads: Number of total reads (including dupliates)
        umis: Number of unique reads / UMIs observ ed

    Returns:
        Estimated library size (total number of unique molecules)
    """
    if umis > reads:
        raise ValueError()
    if umis <= 0:
        raise ValueError()
    # m and M are the lower and upper bound on the saturation, i.e.
    # the fraction of total unique molecules already observed (umis)
    m = 1.0
    M = 100.0
    # Set upper bound
    while coverage_equation(M * umis, umis, reads) > 0:
        M *= 10
    # Find root of coverage_equation
    for _ in range(40):
        r = (m + M) / 2
        u = coverage_equation(r * umis, umis, reads)
        if u > 0:
            m = r
        elif u < 0:
            M = r
        else:
            break
    return umis * (m + M) / 2


def extrapolate_unique(reads: int, umis: int, targetReads: int) -> float:
    """
    Estimate number of unique molecules observed at a given sequencing depth (total reads)
    This is for one 'library', which in our case are the reads from one cell

    Args:
        reads: Genic reads per cell (including duplicates)
        umis: Unique reads (transcript counts)

    Returns:
        Estimated unique reads (counts) at @targetReads (or NaN)
    """
    # Estimate the total library size
    estimatedLibrarySize = estLibSize(reads, umis)
    # Estimate the number of unique reads expected when sequencing targetReads from estLibSize unique molecules
    # (Sampling with replacement)
    res = estimatedLibrarySize * (1 - (((estimatedLibrarySize - 1) / (estimatedLibrarySize)) ** targetReads))
    return int(res) if not np.isnan(res) else res


def calculate_weighted_average(
    weights: list,
    values: list,
    precision: int,
    msg="No passing reads. Please check your sample table or input runfolder / fastq files",
) -> float:
    """
    Calculate the weighted average of a list of values

    Args:
        weights: List of weights for each sample
        values: List of values to average
        precision: Number of decimal places to round to

    Returns:
        Weighted average of the values
    """
    if sum(weights) == 0:
        raise ZeroDivisionError(msg)
    return round(sum(weights[i] * values[i] for i in range(len(weights))) / sum(weights), precision)
