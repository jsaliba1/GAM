"""Summarize mutation-level performance across all LMM output files."""

from __future__ import annotations

import argparse
from pathlib import Path

from gam.logging_utils import configure_logging
from lmm.summarize_lmm_results_folder import summarize_lmm_results_folder


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Summarize mutation calls across an LMM output folder"
    )
    parser.add_argument("--truth-file", type=Path, default=Path("SNP_correct.csv"))
    parser.add_argument("--lmm-folder", type=Path, default=Path("LMM"))
    parser.add_argument("--output", type=Path, default=Path("snpcount_summary.csv"))
    parser.add_argument("--log-level", default="INFO")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    configure_logging(args.log_level)
    summary = summarize_lmm_results_folder(
        truth_file=args.truth_file,
        lmm_folder=args.lmm_folder,
        mode="mutation",
    )
    summary.to_csv(args.output, index=False)
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
