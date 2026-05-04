"""Run TB-specific fine analysis workflow."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from fine_analysis_model_tb import FineAnalysisModel

from gam.logging_utils import configure_logging


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run TB FineAnalysis on curated resistance genes")
    parser.add_argument("--mutation-file", type=Path, default=Path("MUTATIONS.csv"))
    parser.add_argument("--drug-file", type=Path, default=Path("ENA_drug.csv"))
    parser.add_argument("--input-file", type=Path, default=Path("input_max.csv"))
    parser.add_argument("--reg-folder", type=Path, default=Path("Reg"))
    parser.add_argument("--run-num", type=int, default=1)
    parser.add_argument("--output", type=Path, default=Path("fine_analysis_summary.csv"))
    parser.add_argument("--log-level", default="INFO")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    configure_logging(args.log_level)

    analysis = FineAnalysisModel.from_csv(
        mutation_file=args.mutation_file,
        drug_file=args.drug_file,
        input_file=args.input_file,
        reg_folder=args.reg_folder,
        run_num=args.run_num,
    )
    analysis.analyze_snps()
    summary = analysis.score_snps()
    pd.DataFrame(summary).T.to_csv(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
