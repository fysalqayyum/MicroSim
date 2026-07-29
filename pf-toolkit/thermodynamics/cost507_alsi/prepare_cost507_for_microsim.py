#!/usr/bin/env python3
"""Prepare pycalphad's corrected COST507 TDB for the current Al-Si extractor.

Only phase-name and keyword aliases are changed. No FUNCTION or PARAMETER
command, numerical coefficient, temperature range, or constituent model is
altered.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


EXPECTED_UPSTREAM_SHA256 = (
    "6565f4d67695fbbf779d68b9eeeb83e86e44a1a3dde9bafd612d99daffb9c34d"
)
EXPECTED_OUTPUT_SHA256 = (
    "bfba3f8bae2b8f086792cb2ec758b70eff0ff5c2064bfe5394874f489d8a738d"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_tdb", type=Path)
    parser.add_argument("output_tdb", type=Path)
    args = parser.parse_args()

    input_digest = sha256(args.input_tdb)
    if input_digest != EXPECTED_UPSTREAM_SHA256:
        raise ValueError(
            "Unexpected upstream COST507 file: "
            f"expected {EXPECTED_UPSTREAM_SHA256}, found {input_digest}"
        )

    text = args.input_tdb.read_text(encoding="utf-8")
    converted = text.replace("LIQUID:L", "LIQUID")
    converted = converted.replace("DIAMOND_A4", "SI_DIAMOND_A4")
    converted = "\n".join(
        ("CONSTITUENT " + line[6:]) if line.startswith("CONST ") else line
        for line in converted.splitlines()
    )
    if text.endswith("\n"):
        converted += "\n"

    args.output_tdb.write_text(converted, encoding="utf-8")
    output_digest = sha256(args.output_tdb)
    if output_digest != EXPECTED_OUTPUT_SHA256:
        raise RuntimeError(
            "Compatibility conversion was not byte-reproducible: "
            f"expected {EXPECTED_OUTPUT_SHA256}, found {output_digest}"
        )
    print(f"Wrote {args.output_tdb}")
    print(f"SHA-256: {output_digest}")


if __name__ == "__main__":
    main()
