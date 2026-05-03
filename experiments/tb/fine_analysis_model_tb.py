"""Fine-grained SNP analysis model with typed interfaces."""

from __future__ import annotations

import collections
import logging
import math
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd
import scipy.stats as stats

LOGGER = logging.getLogger(__name__)

TARGET_GENES = ["rrs", "rpoB", "katG", "embB", "fabG1", "gyrA"]


def _safe_count(values: pd.Series, key: int) -> int:
    return int(values[values.index == key].iloc[0]) if key in values.index else 0


def _safe_divide(numerator: float, denominator: float) -> float:
    if denominator == 0:
        return float("nan")
    return numerator / denominator


def _binomial_ci(value: float, n: float) -> tuple[float, float]:
    if n == 0 or np.isnan(value):
        return float("nan"), float("nan")
    se = 1.96 * math.sqrt((value * (1 - value)) / n)
    return value - se, value + se


class FineAnalysisModel:
    """Evaluate SNP-level diagnostic performance for a curated gene set."""

    REQUIRED_MUT_COLUMNS = {"UNIQUEID", "MUTATION", "GENE"}
    REQUIRED_INPUT_COLUMNS = {"ID"}

    def __init__(
        self,
        mutations: pd.DataFrame,
        drug_profiles: pd.DataFrame,
        input_isolates: pd.DataFrame,
        reg_folder: Path,
        run_num: int,
    ) -> None:
        missing_mut_cols = self.REQUIRED_MUT_COLUMNS - set(mutations.columns)
        missing_input_cols = self.REQUIRED_INPUT_COLUMNS - set(input_isolates.columns)
        if missing_mut_cols:
            raise ValueError(f"Mutation file is missing columns: {sorted(missing_mut_cols)}")
        if missing_input_cols:
            raise ValueError(f"Input file is missing columns: {sorted(missing_input_cols)}")

        self.mutations = mutations.copy()
        self.drug_profiles = drug_profiles.copy()
        self.folder_names = list(input_isolates["ID"])
        self.reg_folder = reg_folder
        self.run_num = run_num

        self.snp_matrix: Optional[pd.DataFrame] = None
        self.summary: Dict[str, Dict[str, float]] = {}

    @classmethod
    def from_csv(
        cls,
        mutation_file: Path,
        drug_file: Path,
        input_file: Path,
        reg_folder: Path,
        run_num: int = 0,
    ) -> "FineAnalysisModel":
        mutations = pd.read_csv(mutation_file)
        drug_profiles = pd.read_csv(drug_file, index_col="UNIQUEID").replace(["S", "R", "I"], [1, -1, 0])
        input_isolates = pd.read_csv(input_file)
        return cls(
            mutations=mutations,
            drug_profiles=drug_profiles,
            input_isolates=input_isolates,
            reg_folder=reg_folder,
            run_num=run_num,
        )

    def analyze_snps(self) -> pd.DataFrame:
        """Build isolate x mutation representation used for downstream scoring."""
        isolate_snps: dict[str, pd.Series] = {}

        for folder_name in self.folder_names:
            snps = self.mutations[self.mutations["UNIQUEID"] == folder_name]
            snps = snps[snps["GENE"].isin(TARGET_GENES)]
            grouped = snps.groupby(["MUTATION", "GENE"])["MUTATION"].count()
            isolate_snps[folder_name] = grouped

            if self.run_num != 0 and not snps.empty:
                counts = collections.Counter(list(snps["GENE"]))
                repeat_genes = [gene for gene, count in counts.items() if count > 1]
                if counts.get("gyrA", 0) <= 2 and "gyrA" in repeat_genes:
                    repeat_genes.remove("gyrA")

                if repeat_genes:
                    duplicate_snps = snps[snps["GENE"].isin(repeat_genes)]
                    duplicate_idx = duplicate_snps.groupby(["MUTATION", "GENE"])["MUTATION"].count().index
                    for idx in duplicate_idx:
                        isolate_snps[folder_name].at[idx] = 0

        snp_matrix = pd.DataFrame(isolate_snps).replace([1, np.nan], [2, 1])
        snp_matrix = snp_matrix.reindex(sorted(snp_matrix.columns), axis=1)
        self.snp_matrix = snp_matrix
        LOGGER.info("Prepared SNP matrix with shape %s", snp_matrix.shape)
        return snp_matrix

    def score_snps(self) -> Dict[str, Dict[str, float]]:
        """Score SNP associations by drug and summarize performance metrics."""
        if self.snp_matrix is None:
            raise RuntimeError("analyze_snps must be called before score_snps")

        drug_names = [
            "AMIKACIN",
            "ETHAMBUTOL",
            "ETHIONAMIDE",
            "ISONIAZID",
            "KANAMYCIN",
            "LEVOFLOXACIN",
            "MOXIFLOXACIN",
            "RIFAMPICIN",
        ]
        summary: Dict[str, Dict[str, float]] = {}

        for drug_name in drug_names:
            if drug_name not in self.drug_profiles.columns:
                LOGGER.warning("Drug column %s not found; skipping", drug_name)
                continue

            z = self.snp_matrix.mul(self.drug_profiles[drug_name], axis=1)
            reg_file = self.reg_folder / f"{drug_name}.csv"
            reg = pd.read_csv(reg_file, index_col=["MUTATION", "GENE"])
            z = z[z.index.isin(reg.index)]

            meta_rows: Dict[Any, Dict[str, float]] = {}
            for row_index, row in z.iterrows():
                counts = row.value_counts()
                a = _safe_count(counts, -2)
                b = _safe_count(counts, 2)
                c = _safe_count(counts, -1)
                d = _safe_count(counts, 1)

                odds_ratio, p_value = stats.fisher_exact([[a, b], [c, d]])
                sensitivity = _safe_divide(a, a + c)
                specificity = _safe_divide(d, d + b)
                ppv = _safe_divide(a, a + b)
                lci, uci = _binomial_ci(ppv, a + b)

                meta_rows[row_index] = {
                    "a": float(a),
                    "b": float(b),
                    "c": float(c),
                    "d": float(d),
                    "sensitivity": sensitivity,
                    "specificity": specificity,
                    "PPV": ppv,
                    "LCI": lci,
                    "UCI": uci,
                    "OD": odds_ratio,
                    "pval": p_value,
                    "n1": float(a + b),
                }

            meta = pd.DataFrame(meta_rows).T
            if meta.empty:
                LOGGER.warning("No SNP rows available for %s", drug_name)
                continue

            meta["Q"] = meta["n1"] / float(meta["n1"].sum()) * 100
            meta = meta[meta["n1"] >= 5]
            meta = meta[meta["LCI"] >= 0.25]
            meta = meta[meta["OD"] >= 1]

            if len(meta) > 0:
                bonferroni_threshold = 0.05 / len(meta)
                meta = meta[meta["pval"] <= bonferroni_threshold]

            v = reg[reg.index.isin(meta.index)]
            a = b = c = d = e = 0
            for column in v.columns:
                total = float(v[column].sum())
                if total == len(v[column]):
                    d += 1
                elif total == -len(v[column]):
                    c += 1
                elif total > len(v[column]):
                    b += 1
                elif total < -len(v[column]):
                    a += 1
                else:
                    e += 1

            sensitivity = _safe_divide(a, a + c)
            sen_lci, sen_uci = _binomial_ci(sensitivity, a + c)
            specificity = _safe_divide(d, b + d)
            spec_lci, spec_uci = _binomial_ci(specificity, b + d)
            ppv = _safe_divide(a, a + b)
            lci, uci = _binomial_ci(ppv, a + b)

            summary[drug_name] = {
                "a": float(a),
                "b": float(b),
                "c": float(c),
                "d": float(d),
                "sensitivity": sensitivity,
                "sen_LCI": sen_lci,
                "sen_UCI": sen_uci,
                "specificity": specificity,
                "spec_LCI": spec_lci,
                "spec_UCI": spec_uci,
                "PPV": ppv,
                "LCI": lci,
                "UCI": uci,
                "ambiguous_columns": float(e),
            }

        self.summary = summary
        LOGGER.info("Finished FineAnalysis scoring for %s drugs", len(summary))
        return summary
