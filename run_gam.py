"""Run the generic GAM pipeline from a YAML config."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import List, Optional

PROJECT_ROOT = Path(__file__).resolve().parent
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from gam.config import build_gam_run_config, load_yaml
from gam.logging_utils import configure_logging
from gam.models.group_association import GroupAssociationModel


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run generic Group Association Model analysis")
    parser.add_argument("--config", type=Path, required=True, help="YAML config path")
    parser.add_argument("--log-level", default="INFO")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    configure_logging(args.log_level)

    config_dict = load_yaml(args.config)
    config = build_gam_run_config(config_dict)

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
    scores = analysis.score_resistance()
    scores.to_csv(config.output_file)
    analysis.export_profile_pvalues(config.profile_output_file)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
