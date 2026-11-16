"""Diversity and similarity measures used in the paper."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .core import DETECTION_THRESHOLD


def inverse_simpson(profile: pd.Series, threshold: float = DETECTION_THRESHOLD) -> float:
    """Effective number of barcodes, 1 / sum(p_i^2), over detected barcodes."""
    p = profile[profile > threshold].values
    if len(p) == 0:
        return np.nan
    p = p / p.sum()
    return float(1.0 / np.sum(p ** 2))


def cosine_similarity(a: pd.Series, b: pd.Series,
                      threshold: float = DETECTION_THRESHOLD) -> float:
    """Cosine similarity of two barcode frequency vectors.

    Restricted to barcodes above threshold in at least one of the two samples.
    Preferred over Bray-Curtis or Jensen-Shannon because its dependence on the
    number of barcodes in the reference sample is weaker, which makes
    comparisons between regions of differing diversity less confounded.
    """
    joined = pd.concat([a, b], axis=1).fillna(0.0)
    joined.columns = ["a", "b"]
    keep = (joined["a"] > threshold) | (joined["b"] > threshold)
    x, y = joined.loc[keep, "a"].values, joined.loc[keep, "b"].values
    denominator = np.linalg.norm(x) * np.linalg.norm(y)
    return float(x @ y / denominator) if denominator > 0 else np.nan


def expected_barcodes(library: pd.Series, n_cells: int) -> tuple[float, float]:
    """Expected number of distinct barcodes in a sample of n cells.

    E[K] = sum_b [1 - (1 - f_b)^N] over the empirical frequency distribution of
    the library, with the standard deviation from the Poisson binomial variance.
    Equal frequencies are not assumed.
    """
    f = library.values
    absent = np.exp(n_cells * np.log1p(-f))          # (1 - f)^N, stable for small f
    expected = float(np.sum(1.0 - absent))
    variance = float(np.sum(absent * (1.0 - absent)))
    return expected, float(np.sqrt(variance))


def richness_ratio(dataset, individual: str, upstream: str, downstream: str) -> float:
    """Ratio of barcode numbers between two successive samples."""
    a = dataset.barcode_count(individual, upstream)
    b = dataset.barcode_count(individual, downstream)
    if np.isnan(a) or np.isnan(b) or a == 0:
        return np.nan
    return float(b / a)
