"""Convert amino acid mutation notation from three-letter to one-letter format."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import List, Optional

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from gam.logging_utils import configure_logging
from gam.converters.convert_three_letter_mutation import convert_three_letter_mutation


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Convert mutation style for text values or CSV columns")
    parser.add_argument("--name", nargs="+", default=None, help="Inline mutation names")
    parser.add_argument("--input-file", type=Path, default=None, help="CSV file to convert")
    parser.add_argument("--input-column", default="name")
    parser.add_argument("--output-column", default="converted_name")
    parser.add_argument("--output", type=Path, default=Path("converted_mutations.csv"))
    parser.add_argument("--log-level", default="INFO")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    configure_logging(args.log_level)

    if args.name:
        for entry in args.name:
            print(convert_three_letter_mutation(entry))
        return 0

    if args.input_file is None:
        raise ValueError("Provide --name or --input-file")

    frame = pd.read_csv(args.input_file)
    frame[args.output_column] = frame[args.input_column].map(convert_three_letter_mutation)
    frame.to_csv(args.output, index=False)
    print(frame[[args.input_column, args.output_column]].head())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
