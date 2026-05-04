"""Gene-level LMM summary metrics."""

from __future__ import annotations

import pandas as pd


def summarize_lmm_gene_hits(
    truth: pd.DataFrame,
    lmm_results: pd.DataFrame,
    truth_drug_col: str = "DRUG",
    truth_gene_col: str = "GENE",
) -> pd.DataFrame:
    """Summarize TP/FP gene calls per drug."""
    truth_std = truth.rename(columns={truth_drug_col: "Drug", truth_gene_col: "Gene"}).copy()
    lmm_std = lmm_results.rename(columns={"Pheno": "Drug", "SNP": "Mutation"}).copy()
    lmm_std["Gene"] = lmm_std["Mutation"].str.split("_").str[0]

    unique_drugs = truth_std["Drug"].unique()
    rows: list[dict[str, object]] = []

    for drug in unique_drugs:
        true_genes = set(truth_std[truth_std["Drug"] == drug]["Gene"])
        predicted_genes = set(lmm_std[lmm_std["Drug"] == str(drug).upper()]["Gene"])

        true_positives = true_genes & predicted_genes
        false_positives = predicted_genes - true_genes

        rows.append(
            {
                "Drug": drug,
                "True Positives": len(true_positives),
                "False Positives": len(false_positives),
            }
        )

    return pd.DataFrame(rows)
