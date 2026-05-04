"""Folder-level LMM result aggregation."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import pandas as pd

from lmm.filter_significant import filter_significant
from lmm.summarize_lmm_gene_hits import summarize_lmm_gene_hits
from lmm.summarize_lmm_mutation_hits import summarize_lmm_mutation_hits


def summarize_lmm_results_folder(
    truth_file: Path,
    lmm_folder: Path,
    mode: Literal["mutation", "gene"] = "mutation",
) -> pd.DataFrame:
    """Aggregate LMM summary metrics for every CSV result file in a folder."""
    truth = pd.read_csv(truth_file)
    summary_frames: list[pd.DataFrame] = []

    for lmm_file in sorted(lmm_folder.glob("*.csv")):
        lmm_results = pd.read_csv(lmm_file)
        lmm_filtered = filter_significant(lmm_results=lmm_results)

        if mode == "mutation":
            result = summarize_lmm_mutation_hits(truth=truth, lmm_results=lmm_filtered)
        else:
            result = summarize_lmm_gene_hits(truth=truth, lmm_results=lmm_filtered)

        result.insert(0, "File", lmm_file.name)
        summary_frames.append(result)

    if not summary_frames:
        return pd.DataFrame()
    return pd.concat(summary_frames, ignore_index=True)
