"""Run LMM using the TB experiment YAML configuration."""

from __future__ import annotations

import argparse
from pathlib import Path

from gam.logging_utils import configure_logging
from lmm.config import build_lmm_run_config, load_yaml
from lmm.run_lmm_pipeline import run_lmm_pipeline

DEFAULT_CONFIG = Path("experiments/tb/configs/lmm_tb.yaml")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run TB LMM experiment")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--log-level", default="INFO")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    configure_logging(args.log_level)
    config = build_lmm_run_config(load_yaml(args.config))
    run_lmm_pipeline(config=config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
