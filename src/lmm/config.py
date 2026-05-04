"""Configuration models and YAML parsing for generic LMM runs."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class LmmSchema:
    """Column mapping used by the generic LMM pipeline."""

    mutation_sample_id_column: str = "UNIQUEID"
    mutation_gene_column: str = "GENE"
    mutation_name_column: str = "MUTATION"
    phenotype_sample_id_column: str = "UNIQUEID"
    isolate_sample_id_column: str = "ID"
    metadata_sample_id_column: str = "UNIQUEID"
    phenotype_columns: list[str] = field(default_factory=list)
    covariate_columns: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class LmmRunConfig:
    """Run-time settings for generic LMM execution."""

    mutation_file: Path
    phenotype_file: Path
    output_dir: Path
    isolate_file: Path | None = None
    metadata_file: Path | None = None
    sample_sizes: list[int] = field(default_factory=list)
    seed_count: int = 1
    phenotype_value_map: dict[str, int] = field(default_factory=dict)
    schema: LmmSchema = field(default_factory=LmmSchema)


def load_yaml(path: Path) -> dict[str, Any]:
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


def build_lmm_run_config(data: dict[str, Any]) -> LmmRunConfig:
    """Build a typed LMM run config from YAML data."""
    schema_data = data.get("schema", {}) or {}
    run_data = data.get("run", {}) or {}
    input_data = data.get("inputs", {}) or {}
    output_data = data.get("outputs", {}) or {}

    schema = LmmSchema(
        mutation_sample_id_column=schema_data.get("mutation_sample_id_column", "UNIQUEID"),
        mutation_gene_column=schema_data.get("mutation_gene_column", "GENE"),
        mutation_name_column=schema_data.get("mutation_name_column", "MUTATION"),
        phenotype_sample_id_column=schema_data.get("phenotype_sample_id_column", "UNIQUEID"),
        isolate_sample_id_column=schema_data.get("isolate_sample_id_column", "ID"),
        metadata_sample_id_column=schema_data.get("metadata_sample_id_column", "UNIQUEID"),
        phenotype_columns=[str(item) for item in schema_data.get("phenotype_columns", [])],
        covariate_columns=[str(item) for item in schema_data.get("covariate_columns", [])],
    )

    sample_sizes = [int(value) for value in run_data.get("sample_sizes", [])]
    seed_count = int(run_data.get("seed_count", 1))

    phenotype_value_map_data = run_data.get("phenotype_value_map", {"S": 0, "R": 1})
    phenotype_value_map = {str(key): int(value) for key, value in phenotype_value_map_data.items()}

    isolate_file_value = input_data.get("isolate_file")
    metadata_file_value = input_data.get("metadata_file")

    return LmmRunConfig(
        mutation_file=Path(input_data["mutation_file"]),
        phenotype_file=Path(input_data["phenotype_file"]),
        output_dir=Path(output_data.get("output_dir", "lmm_output")),
        isolate_file=Path(isolate_file_value) if isolate_file_value else None,
        metadata_file=Path(metadata_file_value) if metadata_file_value else None,
        sample_sizes=sample_sizes,
        seed_count=seed_count,
        phenotype_value_map=phenotype_value_map,
        schema=schema,
    )
