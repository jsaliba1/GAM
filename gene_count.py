"""Summarize gene-level performance across all LMM output files."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import List, Optional

PROJECT_ROOT = Path(__file__).resolve().parent
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from gam.logging_utils import configure_logging
from lmm.summarize_lmm_results_folder import summarize_lmm_results_folder


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Summarize gene calls across an LMM output folder")
    parser.add_argument("--truth-file", type=Path, default=Path("SNP_correct.csv"))
    parser.add_argument("--lmm-folder", type=Path, default=Path("LMM"))
    parser.add_argument("--output", type=Path, default=Path("genecount_summary.csv"))
    parser.add_argument("--log-level", default="INFO")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    configure_logging(args.log_level)
    summary = summarize_lmm_results_folder(
        truth_file=args.truth_file,
        lmm_folder=args.lmm_folder,
        mode="gene",
    )
    summary.to_csv(args.output, index=False)
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
