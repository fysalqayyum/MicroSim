#!/usr/bin/env python3
"""Convert MicroSim legacy VTK frames to validated VTI/PVD and extract metrics."""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy


def step(path: Path) -> int:
    match = re.search(r"_(\d+)\.vtk$", path.name)
    if not match:
        raise ValueError(f"cannot infer step from {path.name}")
    return int(match.group(1))


def get_array(dataset, candidates):
    point_data = dataset.GetPointData()
    names = [point_data.GetArrayName(i) for i in range(point_data.GetNumberOfArrays())]
    for candidate in candidates:
        if candidate in names:
            return candidate, vtk_to_numpy(point_data.GetArray(candidate))
    raise KeyError(f"none of {candidates} found; arrays={names}")


def read_microsim_fragment(path: Path):
    """Read MicroSim's rank-local binary fragment (not a legacy VTK file)."""
    with path.open("rb") as handle:
        try:
            mesh_x = int(handle.readline())
            mesh_y = int(handle.readline())
            mesh_z = int(handle.readline())
        except ValueError as exc:
            raise ValueError(f"invalid MicroSim fragment header in {path}") from exc
    points = mesh_x * mesh_y * mesh_z
    field_names = ["FCC_A1", "LIQUID", "Mu_Si", "Composition_Si", "Temperature"]
    fields = []
    with path.open("rb") as handle:
        handle.readline()
        handle.readline()
        handle.readline()
        for index, name in enumerate(field_names):
            # Rank-local binary output is big-endian. Convert explicitly to
            # native order before handing the buffer to VTK.
            field = np.fromfile(handle, dtype=">f8", count=points).astype(np.float64)
            if field.size != points:
                raise ValueError(
                    f"{path}: incomplete {name} field ({field.size}/{points})"
                )
            fields.append(field)
            # MicroSim inserts one ASCII newline after every field except the
            # final temperature array.
            if index < len(field_names) - 1:
                separator = handle.read(1)
                if separator != b"\n":
                    raise ValueError(
                        f"{path}: missing binary-field separator after {name}"
                    )

    image = vtk.vtkImageData()
    # MicroSim writes y as the fastest index and uses DIMENSIONS Y X Z in its
    # consolidated VTK writer. Preserve that convention for exact data order.
    image.SetDimensions(mesh_y, mesh_x, mesh_z)
    for name, field in zip(field_names, fields):
        array = numpy_to_vtk(
            np.ascontiguousarray(field),
            deep=True,
        )
        array.SetName(name)
        image.GetPointData().AddArray(array)
    image.GetPointData().SetActiveScalars("FCC_A1")
    return image


def read_frame(path: Path):
    with path.open("rb") as handle:
        signature = handle.read(5)
    if signature == b"# vtk":
        reader = vtk.vtkStructuredPointsReader()
        reader.SetFileName(str(path))
        reader.ReadAllScalarsOn()
        reader.Update()
        return reader.GetOutput()
    return read_microsim_fragment(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir", type=Path)
    parser.add_argument("--prefix", default="")
    parser.add_argument("--dt", type=float, default=5.0e-6)
    parser.add_argument("--dx", type=float, default=7.5e-7)
    parser.add_argument("--shift-file", type=Path)
    parser.add_argument("--series-name", default="AlSi55_G10_V90")
    args = parser.parse_args()

    files = sorted(args.data_dir.glob(f"{args.prefix}*.vtk"), key=step)
    if not files:
        files = sorted(
            args.data_dir.glob(f"Processor_0/{args.prefix}*.vtk"), key=step
        )
    if not files:
        raise SystemExit(f"no VTK frames found in {args.data_dir}")

    out = args.data_dir.parent / "paraview"
    out.mkdir(exist_ok=True)
    metrics = []
    pvd = ET.Element("VTKFile", type="Collection", version="0.1")
    collection = ET.SubElement(pvd, "Collection")

    shift_by_step = {}
    if args.shift_file and args.shift_file.exists():
        shift_by_step = {
            int(t): int(offset)
            for t, offset in np.atleast_2d(np.loadtxt(args.shift_file))
        }

    cumulative_shift = 0
    for source in files:
        current_step = step(source)
        if current_step in shift_by_step:
            cumulative_shift = shift_by_step[current_step]

        data = read_frame(source)
        nx, ny, nz = data.GetDimensions()
        if nz != 1:
            raise ValueError("this post-processor expects a 2D run")

        phase_name, phase_flat = get_array(data, ["alpha", "FCC_A1"])
        comp_name, comp_flat = get_array(
            data, ["Composition_0", "Composition_Si", "Si", "C1"]
        )
        temperature_name, temperature_flat = get_array(data, ["Temperature"])
        if not (
            np.isfinite(phase_flat).all()
            and np.isfinite(comp_flat).all()
            and np.isfinite(temperature_flat).all()
        ):
            raise ValueError(
                f"non-finite field value detected in {source}: "
                f"{phase_name}={np.count_nonzero(~np.isfinite(phase_flat))}, "
                f"{comp_name}={np.count_nonzero(~np.isfinite(comp_flat))}, "
                f"{temperature_name}={np.count_nonzero(~np.isfinite(temperature_flat))}"
            )

        # MicroSim's 2-D writer emits DIMENSIONS MESH_Y MESH_X 1 while its
        # memory loop is x-major/y-minor. After reshaping, the original growth
        # direction y is therefore the second (column) index.
        phase = phase_flat.reshape(ny, nx)
        composition = comp_flat.reshape(ny, nx)
        occupied = np.argwhere(phase >= 0.5)
        tip_local = int(occupied[:, 1].max()) if occupied.size else -1
        tip_global = tip_local + cumulative_shift

        writer = vtk.vtkXMLImageDataWriter()
        target = out / f"{source.stem}.vti"
        writer.SetFileName(str(target))
        writer.SetInputData(data)
        if writer.Write() != 1:
            raise RuntimeError(f"failed to write {target}")

        check = vtk.vtkXMLImageDataReader()
        check.SetFileName(str(target))
        check.Update()
        if check.GetOutput().GetNumberOfPoints() != data.GetNumberOfPoints():
            raise RuntimeError(f"point-count mismatch after writing {target}")

        ET.SubElement(
            collection,
            "DataSet",
            timestep=f"{current_step * args.dt:.12g}",
            group="",
            part="0",
            file=target.name,
        )
        metrics.append(
            {
                "step": current_step,
                "time_s": current_step * args.dt,
                "phase_field": phase_name,
                "composition_field": comp_name,
                "temperature_field": temperature_name,
                "solid_area_fraction": float(np.mean(phase >= 0.5)),
                "tip_local_um": tip_local * args.dx * 1e6,
                "shift_um": cumulative_shift * args.dx * 1e6,
                "tip_global_um": tip_global * args.dx * 1e6,
                "phase_min": float(np.min(phase_flat)),
                "phase_max": float(np.max(phase_flat)),
                "xSi_min": float(np.nanmin(composition)),
                "xSi_mean": float(np.nanmean(composition)),
                "xSi_max": float(np.nanmax(composition)),
                "temperature_min_K": float(np.min(temperature_flat)),
                "temperature_max_K": float(np.max(temperature_flat)),
            }
        )

    times = np.asarray([row["time_s"] for row in metrics], dtype=float)
    tips = np.asarray([row["tip_global_um"] for row in metrics], dtype=float)
    velocities = np.gradient(tips, times) if len(metrics) > 1 else np.asarray([np.nan])
    for row, velocity in zip(metrics, velocities):
        row["tip_velocity_um_s"] = float(velocity)

    tree = ET.ElementTree(pvd)
    ET.indent(tree, space="  ")
    pvd_path = out / f"{args.series_name}.pvd"
    tree.write(pvd_path, encoding="utf-8", xml_declaration=True)
    with (out / "frame_metrics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=metrics[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(metrics)
    print(f"validated {len(files)} frames; ParaView series: {pvd_path}")


if __name__ == "__main__":
    main()
