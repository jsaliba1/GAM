"""Generic Group Association Model implementation."""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import scipy.stats as stats

from gam.config import GamSchema

LOGGER = logging.getLogger(__name__)

Profile = Tuple[str, ...]


@dataclass(frozen=True)
class GroupingResult:
    grouped_drugs: pd.DataFrame
    profiles: List[Profile]


class GroupAssociationModel:
    """Detect mutation associations across grouped resistance profiles."""

    def __init__(
        self,
        mutations: pd.DataFrame,
        drug_profiles: pd.DataFrame,
        bugs: pd.DataFrame,
        schema: Optional[GamSchema] = None,
    ) -> None:
        self.schema = schema or GamSchema()
        self._validate_input_columns(mutations=mutations, bugs=bugs)

        self.mutations = self._filter_synonymous_mutations(mutations)
        self.drug_profiles = self._normalize_drug_profiles(drug_profiles)
        self.bugs = bugs.copy()
        self.drug_names = self.drug_profiles.columns.tolist()
        self.n_drugs = len(self.drug_names)

        self.grouped_drugs: Optional[pd.DataFrame] = None
        self.profile_order: List[Profile] = []
        self.p_values: Dict[Profile, pd.Series] = {}

    @classmethod
    def from_csv(
        cls,
        mutation_file: Path,
        drug_file: Path,
        bug_file: Path,
        schema: Optional[GamSchema] = None,
    ) -> "GroupAssociationModel":
        """Create a model instance from CSV input paths."""
        schema_to_use = schema or GamSchema()

        mutations = pd.read_csv(mutation_file, low_memory=False)
        drug_profiles = pd.read_csv(drug_file)
        if schema_to_use.drug_sample_id_column not in drug_profiles.columns:
            raise ValueError(
                f"Drug file is missing sample id column: {schema_to_use.drug_sample_id_column}"
            )
        drug_profiles = drug_profiles.set_index(schema_to_use.drug_sample_id_column)

        bugs = pd.read_csv(bug_file)
        return cls(mutations=mutations, drug_profiles=drug_profiles, bugs=bugs, schema=schema_to_use)

    def group_bugs(self, mask_probability: float = 0.0, remove_n: int = 0, seed: int = 0) -> GroupingResult:
        """Mask phenotype profiles and group isolates by profile."""
        if not 0.0 <= mask_probability <= 1.0:
            raise ValueError("mask_probability must be between 0 and 1")
        if remove_n < 0:
            raise ValueError("remove_n cannot be negative")
        if remove_n > len(self.bugs):
            raise ValueError("remove_n cannot exceed number of bug rows")

        rng = np.random.default_rng(seed)
        if remove_n == 0:
            bugs_subset = self.bugs
        else:
            drop_indices = rng.choice(self.bugs.index.to_numpy(), remove_n, replace=False)
            bugs_subset = self.bugs.drop(drop_indices)

        bug_names = bugs_subset[self.schema.bug_sample_id_column].astype(str).tolist()
        drugs = self.drug_profiles[self.drug_profiles.index.isin(bug_names)].copy()

        mask = rng.choice([0, 1], size=drugs.shape, p=[mask_probability, 1 - mask_probability])
        grouped_drugs = drugs.where(mask == 1, other="a")

        profiles = [tuple(str(x) for x in row) for _, row in grouped_drugs.iterrows()]
        counts = pd.Series(Counter(profiles))
        counts = counts.sort_values(ascending=False)

        self.grouped_drugs = grouped_drugs
        self.profile_order = list(counts.index)
        LOGGER.info("Created %s unique resistance profiles", len(self.profile_order))
        return GroupingResult(grouped_drugs=grouped_drugs, profiles=self.profile_order)

    def analyze_resistance(self) -> Dict[Profile, pd.Series]:
        """Run profile-wise Fisher tests versus the all-susceptible control profile."""
        if self.grouped_drugs is None or not self.profile_order:
            raise RuntimeError("group_bugs must be called before analyze_resistance")

        control_profile = tuple(self.schema.susceptible_label for _ in range(self.n_drugs))
        control_ids = self._strain_ids_for_profile(control_profile)
        if len(control_ids) <= 1:
            LOGGER.warning("Control group not found or too small; skipping resistance analysis")
            self.p_values = {}
            return self.p_values

        control_counts = self._mutation_counts(control_ids)
        p_values: Dict[Profile, pd.Series] = {}

        for profile in self.profile_order:
            if profile == control_profile:
                continue

            group_ids = self._strain_ids_for_profile(profile)
            if len(group_ids) <= 1:
                continue

            resistant_counts = self._mutation_counts(group_ids)
            combined = pd.concat([resistant_counts, control_counts], axis=1).fillna(0)
            combined.columns = ["resistant", "control"]
            combined = combined[combined["resistant"] > 0]
            if combined.empty:
                continue

            bonferroni_threshold = 0.05 / len(combined)
            significant: List[Tuple[Tuple[str, str], float]] = []

            for row_index, row in combined.iterrows():
                a = int(row["resistant"])
                b = int(row["control"])
                contingency = [[a, b], [len(group_ids) - a, len(control_ids) - b]]
                odds_ratio, p_value = stats.fisher_exact(contingency)
                if p_value <= bonferroni_threshold and odds_ratio > 1:
                    significant.append((row_index, p_value))

            if significant:
                idx = pd.MultiIndex.from_tuples(
                    [item[0] for item in significant],
                    names=[self.schema.mutation_name_column, self.schema.mutation_gene_column],
                )
                series = pd.Series([item[1] for item in significant], index=idx, dtype=float)
                p_values[profile] = series

        self.p_values = p_values
        LOGGER.info("Found significant mutations in %s profiles", len(p_values))
        return self.p_values

    def score_resistance(self) -> pd.DataFrame:
        """Summarize mutation-level profile significance into drug-level p-values."""
        if not self.p_values:
            return pd.DataFrame(columns=self.drug_names)

        pval_df = pd.DataFrame(self.p_values).fillna(1.0)
        md_pval = pd.DataFrame(1.0, index=pval_df.index, columns=self.drug_names)

        for mutation_idx, mut_pvals in pval_df.iterrows():
            a = [0] * self.n_drugs
            b = [0] * self.n_drugs
            c = [0] * self.n_drugs
            d = [0] * self.n_drugs

            for profile, p_value in mut_pvals.items():
                for idx, resistance_state in enumerate(profile):
                    if p_value < 0.05:
                        if resistance_state == self.schema.susceptible_label:
                            c[idx] += 1
                        elif resistance_state == self.schema.resistant_label:
                            a[idx] += 1
                    else:
                        if resistance_state == self.schema.susceptible_label:
                            d[idx] += 1
                        elif resistance_state == self.schema.resistant_label:
                            b[idx] += 1

            for idx, drug_name in enumerate(self.drug_names):
                contingency = [[a[idx], b[idx]], [c[idx], d[idx]]]
                odds_ratio, p_value = stats.fisher_exact(contingency)
                md_pval.at[mutation_idx, drug_name] = p_value if odds_ratio > 1 else 1.0

        bonferroni_threshold = 0.05 / len(md_pval)
        filtered = md_pval.where(md_pval <= bonferroni_threshold)
        filtered = filtered.dropna(axis=0, how="all")
        return filtered

    def export_profile_pvalues(self, output_file: Path) -> None:
        """Persist profile-level p-values to CSV."""
        if not self.p_values:
            pd.DataFrame().to_csv(output_file)
            return
        output = pd.DataFrame(self.p_values).fillna(1.0)
        output.to_csv(output_file)

    def _filter_synonymous_mutations(self, mutations: pd.DataFrame) -> pd.DataFrame:
        synonymous_col = self.schema.mutation_is_synonymous_column
        if synonymous_col is None or synonymous_col not in mutations.columns:
            return mutations.copy()
        return mutations[mutations[synonymous_col] == False].copy()  # noqa: E712

    def _normalize_drug_profiles(self, drug_profiles: pd.DataFrame) -> pd.DataFrame:
        normalized = drug_profiles.copy()
        normalized.index = normalized.index.astype(str)
        normalized = normalized.fillna("a").astype(str)
        return normalized

    def _validate_input_columns(self, mutations: pd.DataFrame, bugs: pd.DataFrame) -> None:
        required_mut_cols = {
            self.schema.mutation_sample_id_column,
            self.schema.mutation_name_column,
            self.schema.mutation_gene_column,
        }

        missing_mut_cols = required_mut_cols - set(mutations.columns)
        missing_bug_cols = {self.schema.bug_sample_id_column} - set(bugs.columns)

        if missing_mut_cols:
            raise ValueError(f"Mutation file is missing columns: {sorted(missing_mut_cols)}")
        if missing_bug_cols:
            raise ValueError(f"Bug file is missing columns: {sorted(missing_bug_cols)}")

    def _strain_ids_for_profile(self, profile: Profile) -> List[str]:
        if self.grouped_drugs is None:
            raise RuntimeError("group_bugs must be called before querying profile groups")

        filtered = self.grouped_drugs
        for drug_name, status in zip(self.drug_names, profile):
            filtered = filtered[filtered[drug_name] == status]
        return list(filtered.index.astype(str))

    def _mutation_counts(self, strain_ids: Sequence[str]) -> pd.Series:
        subset = self.mutations[self.mutations[self.schema.mutation_sample_id_column].astype(str).isin(strain_ids)]
        if subset.empty:
            return pd.Series(dtype=float)
        grouped = subset.groupby(
            [self.schema.mutation_name_column, self.schema.mutation_gene_column]
        )[self.schema.mutation_name_column].count()
        return grouped


def format_profile(profile: Iterable[str]) -> str:
    """Build a human-readable profile label."""
    return "[" + ", ".join(profile) + "]"
