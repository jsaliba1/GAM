"""Convert amino acid substitutions from three-letter to one-letter notation."""

from __future__ import annotations

import re

THREE_TO_ONE = {
    "VAL": "V",
    "ILE": "I",
    "LEU": "L",
    "GLU": "E",
    "GLN": "Q",
    "ASP": "D",
    "ASN": "N",
    "HIS": "H",
    "TRP": "W",
    "PHE": "F",
    "TYR": "Y",
    "ARG": "R",
    "LYS": "K",
    "SER": "S",
    "THR": "T",
    "MET": "M",
    "ALA": "A",
    "GLY": "G",
    "PRO": "P",
    "CYS": "C",
}


def convert_three_letter_mutation(name: str) -> str:
    """Convert names like `Ser450Leu` into one-letter form (`S450L`)."""
    parts = re.split(r"(\d+)", name)
    if len(parts) < 3:
        return name

    left, middle, right = parts[0], parts[1], parts[2]
    left_up = left.upper()
    right_up = right.upper()

    if len(left_up) == 3 and len(right_up) == 3 and right_up not in {"DEL", "DUP"}:
        return f"{THREE_TO_ONE[left_up]}{middle}{THREE_TO_ONE[right_up]}"
    if len(left_up) == 3 and right == "*":
        return f"{THREE_TO_ONE[left_up]}{middle}!"
    if len(left_up) == 3:
        return f"{THREE_TO_ONE[left_up]}{middle}{right}"
    return f"{left}{middle}{right}"
