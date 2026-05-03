from gam.config import build_gam_run_config
from lmm.config import build_lmm_run_config


def test_build_gam_run_config_uses_custom_schema() -> None:
    config = build_gam_run_config(
        {
            "inputs": {
                "mutation_file": "mut.csv",
                "drug_file": "drug.csv",
                "bug_file": "bugs.csv",
            },
            "outputs": {"scores_file": "scores.csv", "profile_pvalues_file": "profiles.csv"},
            "schema": {
                "mutation_sample_id_column": "sample",
                "bug_sample_id_column": "isolate",
                "susceptible_label": "0",
                "resistant_label": "1",
            },
        }
    )
    assert str(config.mutation_file) == "mut.csv"
    assert config.schema.mutation_sample_id_column == "sample"
    assert config.schema.bug_sample_id_column == "isolate"
    assert config.schema.susceptible_label == "0"


def test_build_lmm_run_config_reads_covariates() -> None:
    config = build_lmm_run_config(
        {
            "inputs": {
                "mutation_file": "mut.csv",
                "phenotype_file": "pheno.csv",
                "isolate_file": "isolates.csv",
                "metadata_file": "meta.csv",
            },
            "outputs": {"output_dir": "out/lmm"},
            "run": {"sample_sizes": [50], "seed_count": 3, "phenotype_value_map": {"S": 0, "R": 1}},
            "schema": {"covariate_columns": ["lineage", "site"]},
        }
    )
    assert str(config.output_dir) == "out/lmm"
    assert config.sample_sizes == [50]
    assert config.seed_count == 3
    assert config.schema.covariate_columns == ["lineage", "site"]
