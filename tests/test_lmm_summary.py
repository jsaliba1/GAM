import pandas as pd

from lmm.filter_significant import filter_significant
from lmm.summarize_lmm_mutation_hits import summarize_lmm_mutation_hits


def test_filter_significant_uses_bonferroni() -> None:
    lmm = pd.DataFrame(
        {
            "PValue": [0.0001, 0.02, 0.03, 0.0002],
            "Pheno": ["A", "A", "B", "B"],
            "SNP": ["g1_a", "g2_b", "g3_c", "g4_d"],
        }
    )
    filtered = filter_significant(lmm_results=lmm, alpha=0.05)
    assert len(filtered) == 2


def test_summarize_lmm_mutation_hits_counts_tp_fp_fn() -> None:
    truth = pd.DataFrame({"DRUG": ["DrugA", "DrugA"], "MUT": ["g1_a", "g2_b"]})
    lmm = pd.DataFrame(
        {
            "Pheno": ["DRUGA", "DRUGA"],
            "SNP": ["g1_a", "g9_z"],
            "PValue": [1e-8, 1e-8],
        }
    )
    summary = summarize_lmm_mutation_hits(truth=truth, lmm_results=lmm)
    row = summary.iloc[0]
    assert int(row["True Positives"]) == 1
    assert int(row["False Positives"]) == 1
    assert int(row["False Negatives"]) == 1
