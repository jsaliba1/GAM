"""Bonferroni filtering for LMM output tables."""

from __future__ import annotations

import pandas as pd


def filter_significant(lmm_results: pd.DataFrame, alpha: float = 0.05) -> pd.DataFrame:
    """Apply Bonferroni correction and keep significant rows only."""
    if lmm_results.empty:
        return lmm_results.copy()
    threshold = alpha / len(lmm_results)
    return lmm_results[lmm_results["PValue"] < threshold].copy()
