"""Run GAM using the TB experiment YAML configuration."""

from __future__ import annotations

import argparse
from pathlib import Path

from gam.config import build_gam_run_config, load_yaml
from gam.logging_utils import configure_logging
from gam.models.group_association import GroupAssociationModel

DEFAULT_CONFIG = Path("experiments/tb/configs/gam_tb.yaml")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run TB GAM experiment")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--log-level", default="INFO")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    configure_logging(args.log_level)

    config = build_gam_run_config(load_yaml(args.config))
    analysis = GroupAssociationModel.from_csv(
        mutation_file=config.mutation_file,
        drug_file=config.drug_file,
        bug_file=config.bug_file,
        schema=config.schema,
    )
    analysis.group_bugs(
        mask_probability=config.mask_probability,
        remove_n=config.remove_n,
        seed=config.seed,
    )
    analysis.analyze_resistance()
    analysis.score_resistance().to_csv(config.output_file)
    analysis.export_profile_pvalues(config.profile_output_file)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
