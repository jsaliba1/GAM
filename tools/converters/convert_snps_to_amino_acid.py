"""Convert nucleotide SNP edits to amino-acid substitutions."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import List, Optional

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from gam.converters.convert_snps_to_aa import convert_snps_to_aa
from gam.logging_utils import configure_logging


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Convert SNP edits to amino-acid notation")
    parser.add_argument("--sequence", required=True, help="Reference DNA sequence")
    parser.add_argument("--strand", choices=["+", "-"], required=True)
    parser.add_argument("--start-pos", type=int, required=True)
    parser.add_argument("--end-pos", type=int, required=True)
    parser.add_argument("--positions", required=True, help="Comma separated genomic positions")
    parser.add_argument("--refs", required=True, help="Comma separated reference bases")
    parser.add_argument("--alts", required=True, help="Comma separated alternative bases")
    parser.add_argument("--log-level", default="INFO")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    configure_logging(args.log_level)

    positions = [int(x.strip()) for x in args.positions.split(",") if x.strip()]
    refs = [x.strip() for x in args.refs.split(",") if x.strip()]
    alts = [x.strip() for x in args.alts.split(",") if x.strip()]

    mutations = convert_snps_to_aa(
        sequence=args.sequence,
        strand=args.strand,
        start_pos=args.start_pos,
        end_pos=args.end_pos,
        positions=positions,
        refs=refs,
        alts=alts,
    )
    for mutation in mutations:
        print(mutation)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
