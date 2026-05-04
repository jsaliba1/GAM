"""Generic LMM association pipeline."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from lmm.config import LmmRunConfig

LOGGER = logging.getLogger(__name__)


def _make_snp_data(data: pd.DataFrame, snp_data_cls: type) -> object:
    iid = [["fam0", sample_id] for sample_id in data.index]
    sid = data.columns.tolist()
    val = data.values
    return snp_data_cls(iid=iid, sid=sid, val=val)


def _resolve_sample_ids(
    config: LmmRunConfig, mutations: pd.DataFrame, phenotypes: pd.DataFrame
) -> list[str]:
    mutation_ids = set(mutations[config.schema.mutation_sample_id_column].astype(str))
    phenotype_ids = set(phenotypes.index.astype(str))
    shared = sorted(mutation_ids & phenotype_ids)

    if config.isolate_file is None:
        return shared

    isolate_df = pd.read_csv(config.isolate_file)
    requested_ids = isolate_df[config.schema.isolate_sample_id_column].astype(str).tolist()
    return [sample_id for sample_id in requested_ids if sample_id in shared]


def _build_mutation_matrix(
    config: LmmRunConfig, mutations: pd.DataFrame, sample_ids: Sequence[str]
) -> pd.DataFrame:
    filtered = mutations[
        mutations[config.schema.mutation_sample_id_column].astype(str).isin(sample_ids)
    ].copy()
    filtered["composite_snp"] = (
        filtered[config.schema.mutation_gene_column].astype(str)
        + "_"
        + filtered[config.schema.mutation_name_column].astype(str)
    )
    filtered["value"] = 1

    mutation_matrix = filtered.pivot_table(
        index="composite_snp",
        columns=config.schema.mutation_sample_id_column,
        values="value",
        aggfunc="max",
    )
    mutation_matrix = mutation_matrix.fillna(0).astype("int8")
    mutation_matrix = mutation_matrix.reindex(columns=list(sample_ids), fill_value=0)
    return mutation_matrix


def _build_kinship_matrix(mutation_matrix: pd.DataFrame) -> pd.DataFrame:
    mutation_array = mutation_matrix.values
    intersection_matrix = np.dot(mutation_array.T, mutation_array)
    row_sums = mutation_array.sum(axis=0)
    union_matrix = np.add.outer(row_sums, row_sums) - intersection_matrix
    union_matrix[union_matrix == 0] = 1
    kinship_matrix = intersection_matrix / union_matrix
    np.fill_diagonal(kinship_matrix, 1)
    return pd.DataFrame(
        kinship_matrix, index=mutation_matrix.columns, columns=mutation_matrix.columns
    )


def _prepare_phenotypes(
    config: LmmRunConfig, phenotypes: pd.DataFrame, sample_ids: Sequence[str]
) -> pd.DataFrame:
    subset = phenotypes.loc[list(sample_ids)].copy()
    if config.schema.phenotype_columns:
        subset = subset[config.schema.phenotype_columns]

    if config.phenotype_value_map:
        subset = subset.replace(config.phenotype_value_map)

    for column in subset.columns:
        subset[column] = pd.to_numeric(subset[column], errors="coerce")
    subset = subset.dropna(axis=1, how="all")
    if subset.isna().any().any():
        raise ValueError("Phenotype table has non-numeric or missing values after mapping")
    return subset.astype("int8")


def _prepare_covariates(
    config: LmmRunConfig,
    sample_ids: Sequence[str],
) -> pd.DataFrame | None:
    if config.metadata_file is None or not config.schema.covariate_columns:
        return None

    metadata = pd.read_csv(config.metadata_file)
    metadata = metadata.set_index(config.schema.metadata_sample_id_column)
    metadata.index = metadata.index.astype(str)
    metadata = metadata.reindex(list(sample_ids))

    covariate_frames: list[pd.DataFrame] = []
    for column in config.schema.covariate_columns:
        if column not in metadata.columns:
            continue
        dummies = pd.get_dummies(metadata[column].fillna("Unknown"), drop_first=True)
        covariate_frames.append(dummies)

    if not covariate_frames:
        return None
    return pd.concat(covariate_frames, axis=1)


def _resolve_sample_sizes(config: LmmRunConfig, n_samples: int) -> list[int]:
    if config.sample_sizes:
        sizes = config.sample_sizes
    else:
        sizes = [n_samples]
    for size in sizes:
        if size > n_samples:
            raise ValueError(f"Requested sample size {size} exceeds available isolates {n_samples}")
    return sizes


def run_lmm_pipeline(config: LmmRunConfig) -> None:
    """Run LMM association experiments using generic column mappings."""
    from fastlmm.association import single_snp
    from pysnptools.snpreader import SnpData
    from sklearn.utils import resample

    mutations = pd.read_csv(config.mutation_file, low_memory=False)
    phenotypes = pd.read_csv(config.phenotype_file)
    phenotypes = phenotypes.set_index(config.schema.phenotype_sample_id_column)
    phenotypes.index = phenotypes.index.astype(str)

    sample_ids = _resolve_sample_ids(config=config, mutations=mutations, phenotypes=phenotypes)
    if not sample_ids:
        raise ValueError("No shared sample IDs found across mutation and phenotype inputs")

    mutation_matrix = _build_mutation_matrix(
        config=config, mutations=mutations, sample_ids=sample_ids
    )
    kinship_matrix = _build_kinship_matrix(mutation_matrix=mutation_matrix)
    phenotypes_clean = _prepare_phenotypes(
        config=config, phenotypes=phenotypes, sample_ids=sample_ids
    )
    covariates = _prepare_covariates(config=config, sample_ids=sample_ids)

    config.output_dir.mkdir(parents=True, exist_ok=True)
    sample_sizes = _resolve_sample_sizes(config=config, n_samples=len(sample_ids))

    for size in sample_sizes:
        for seed in range(config.seed_count):
            LOGGER.info("Running LMM for sample_size=%s seed=%s", size, seed)
            sampled_ids = resample(sample_ids, n_samples=size, replace=False, random_state=seed)

            phenotype_snp = _make_snp_data(phenotypes_clean.loc[sampled_ids], SnpData)
            snp_matrix = _make_snp_data(mutation_matrix[sampled_ids].T, SnpData)
            kinship_snp = _make_snp_data(kinship_matrix[sampled_ids].loc[sampled_ids], SnpData)
            covariate_snp = None
            if covariates is not None and not covariates.empty:
                covariate_snp = _make_snp_data(covariates.loc[sampled_ids], SnpData)

            kwargs = {
                "test_snps": snp_matrix,
                "pheno": phenotype_snp,
                "K0": kinship_snp,
                "count_A1": False,
                "leave_out_one_chrom": False,
            }
            if covariate_snp is not None:
                kwargs["covar"] = covariate_snp

            results_df = single_snp(**kwargs)
            output_file = Path(config.output_dir) / f"{size}_lmm_{seed}.csv"
            results_df.to_csv(output_file)
