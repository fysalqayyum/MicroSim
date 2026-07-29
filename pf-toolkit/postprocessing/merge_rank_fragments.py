#!/usr/bin/env python3
"""Merge x-decomposed MicroSim rank fragments into one binary fragment."""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

import numpy as np


FIELDS = ("FCC_A1", "LIQUID", "Mu_Si", "Composition_Si", "Temperature")


def rank(path: Path) -> int:
    match = re.fullmatch(r"Processor_(\d+)", path.parent.name)
    if match is None:
        raise ValueError(f"cannot infer rank from {path}")
    return int(match.group(1))


def step(path: Path) -> int:
    match = re.search(r"_(\d+)\.vtk$", path.name)
    if match is None:
        raise ValueError(f"cannot infer step from {path}")
    return int(match.group(1))


def read_fragment(path: Path) -> tuple[tuple[int, int, int], dict[str, np.ndarray]]:
    with path.open("rb") as handle:
        dimensions = tuple(int(handle.readline()) for _ in range(3))
        count = int(np.prod(dimensions))
        arrays: dict[str, np.ndarray] = {}
        for index, field in enumerate(FIELDS):
            values = np.fromfile(handle, dtype=">f8", count=count)
            if values.size != count:
                raise ValueError(f"{path}: incomplete {field}")
            arrays[field] = values.reshape(dimensions)
            if index < len(FIELDS) - 1 and handle.read(1) != b"\n":
                raise ValueError(f"{path}: missing separator after {field}")
    return dimensions, arrays


def write_fragment(path: Path, arrays: dict[str, np.ndarray]) -> None:
    dimensions = arrays[FIELDS[0]].shape
    if any(arrays[field].shape != dimensions for field in FIELDS):
        raise ValueError("merged fields do not have a common shape")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        for dimension in dimensions:
            handle.write(f"{dimension}\n".encode())
        for index, field in enumerate(FIELDS):
            np.asarray(arrays[field], dtype=">f8").tofile(handle)
            if index < len(FIELDS) - 1:
                handle.write(b"\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("data_root", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    grouped: dict[int, list[Path]] = defaultdict(list)
    for path in args.data_root.glob(f"Processor_*/*{args.prefix}*.vtk"):
        grouped[step(path)].append(path)
    if not grouped:
        raise SystemExit(f"no rank fragments found below {args.data_root}")

    for timestep, paths in sorted(grouped.items()):
        paths.sort(key=rank)
        observed_ranks = [rank(path) for path in paths]
        expected_ranks = list(range(len(paths)))
        if observed_ranks != expected_ranks:
            raise ValueError(
                f"step {timestep}: ranks {observed_ranks} != {expected_ranks}"
            )

        pieces = [read_fragment(path) for path in paths]
        yz_shapes = {(dimensions[1], dimensions[2]) for dimensions, _ in pieces}
        if len(yz_shapes) != 1:
            raise ValueError(f"step {timestep}: inconsistent y/z dimensions")
        merged = {
            field: np.concatenate(
                [arrays[field] for _, arrays in pieces], axis=0
            )
            for field in FIELDS
        }
        output = args.output_dir / f"{args.prefix}_{timestep}.vtk"
        write_fragment(output, merged)
        print(
            f"step={timestep} ranks={len(paths)} "
            f"shape={merged[FIELDS[0]].shape} output={output}"
        )


if __name__ == "__main__":
    main()
