"""Configuration models and YAML parsing for generic GAM runs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional


@dataclass(frozen=True)
class GamSchema:
    """Column and label mapping used by the generic GAM pipeline."""

    mutation_sample_id_column: str = "UNIQUEID"
    mutation_name_column: str = "MUTATION"
    mutation_gene_column: str = "GENE"
    mutation_is_synonymous_column: Optional[str] = "IS_SYNONYMOUS"
    bug_sample_id_column: str = "ID"
    drug_sample_id_column: str = "UNIQUEID"
    susceptible_label: str = "S"
    resistant_label: str = "R"


@dataclass(frozen=True)
class GamRunConfig:
    """Run-time settings for generic GAM execution."""

    mutation_file: Path
    drug_file: Path
    bug_file: Path
    output_file: Path
    profile_output_file: Path
    mask_probability: float = 0.0
    remove_n: int = 0
    seed: int = 0
    schema: GamSchema = GamSchema()


def load_yaml(path: Path) -> Dict[str, Any]:
    """Load a YAML dictionary from disk."""
    try:
        import yaml
    except ImportError as exc:  # pragma: no cover - import guard
        raise ImportError("PyYAML is required to read YAML configs") from exc

    with path.open("r", encoding="utf-8") as handle:
        parsed = yaml.safe_load(handle) or {}
    if not isinstance(parsed, dict):
        raise ValueError("YAML config must be a top-level mapping")
    return parsed


def build_gam_run_config(data: Dict[str, Any]) -> GamRunConfig:
    """Build a typed GAM run config from YAML data."""
    schema_data = data.get("schema", {}) or {}
    schema = GamSchema(
        mutation_sample_id_column=schema_data.get("mutation_sample_id_column", "UNIQUEID"),
        mutation_name_column=schema_data.get("mutation_name_column", "MUTATION"),
        mutation_gene_column=schema_data.get("mutation_gene_column", "GENE"),
        mutation_is_synonymous_column=schema_data.get("mutation_is_synonymous_column", "IS_SYNONYMOUS"),
        bug_sample_id_column=schema_data.get("bug_sample_id_column", "ID"),
        drug_sample_id_column=schema_data.get("drug_sample_id_column", "UNIQUEID"),
        susceptible_label=str(schema_data.get("susceptible_label", "S")),
        resistant_label=str(schema_data.get("resistant_label", "R")),
    )

    run_data = data.get("run", {}) or {}
    input_data = data.get("inputs", {}) or {}
    output_data = data.get("outputs", {}) or {}

    return GamRunConfig(
        mutation_file=Path(input_data["mutation_file"]),
        drug_file=Path(input_data["drug_file"]),
        bug_file=Path(input_data["bug_file"]),
        output_file=Path(output_data.get("scores_file", "gam_scores.csv")),
        profile_output_file=Path(output_data.get("profile_pvalues_file", "gam_profile_pvalues.csv")),
        mask_probability=float(run_data.get("mask_probability", 0.0)),
        remove_n=int(run_data.get("remove_n", 0)),
        seed=int(run_data.get("seed", 0)),
        schema=schema,
    )
