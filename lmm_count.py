"""Summarize a single LMM results file against mutation truth labels."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from gam.logging_utils import configure_logging
from lmm.filter_significant import filter_significant
from lmm.summarize_lmm_mutation_hits import summarize_lmm_mutation_hits


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Count true/false mutation calls for one LMM result file"
    )
    parser.add_argument("--truth-file", type=Path, default=Path("SNP_correct.csv"))
    parser.add_argument("--lmm-file", type=Path, default=Path("LMM/7179_lmm_0.csv"))
    parser.add_argument("--output", type=Path, default=Path("lmmcount_summary.csv"))
    parser.add_argument("--log-level", default="INFO")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    configure_logging(args.log_level)

    truth = pd.read_csv(args.truth_file)
    lmm_results = pd.read_csv(args.lmm_file)
    lmm_filtered = filter_significant(lmm_results)
    summary = summarize_lmm_mutation_hits(truth=truth, lmm_results=lmm_filtered)
    summary.to_csv(args.output, index=False)
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
