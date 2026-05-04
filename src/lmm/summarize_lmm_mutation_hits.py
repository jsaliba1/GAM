"""Mutation-level LMM summary metrics."""

from __future__ import annotations

import pandas as pd


def summarize_lmm_mutation_hits(
    truth: pd.DataFrame,
    lmm_results: pd.DataFrame,
    truth_drug_col: str = "DRUG",
    truth_mut_col: str = "MUT",
) -> pd.DataFrame:
    """Summarize TP/FP/FN mutation calls per drug."""
    truth_std = truth.rename(columns={truth_drug_col: "Drug", truth_mut_col: "Mutation"}).copy()
    lmm_std = lmm_results.rename(columns={"Pheno": "Drug", "SNP": "Mutation"}).copy()

    unique_drugs = truth_std["Drug"].unique()
    rows: list[dict[str, object]] = []

    for drug in unique_drugs:
        true_mutations = set(truth_std[truth_std["Drug"] == drug]["Mutation"])
        predicted_mutations = set(lmm_std[lmm_std["Drug"] == str(drug).upper()]["Mutation"])

        true_positives = true_mutations & predicted_mutations
        false_positives = predicted_mutations - true_mutations
        false_negatives = true_mutations - predicted_mutations

        rows.append(
            {
                "Drug": drug,
                "True Positives": len(true_positives),
                "False Positives": len(false_positives),
                "False Negatives": len(false_negatives),
            }
        )

    return pd.DataFrame(rows)
