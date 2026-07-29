#!/usr/bin/env python3
"""
Re-embed a developed MicroSim front into a new domain.

Takes a restart frame from a previous run, crops it in x (with a seam-matched
crop offset), pads liquid above it to a taller domain, and writes a new restart
frame that the solver will accept as step 0.  The thermal field is NOT carried
over -- the solver overwrites every temperature cell from `Tempgrady` in
Input.in, so restarting the output of this script with STARTTIME=0 and
RESTART=1 discards all accumulated cooling.  That discard is the point: it
removes the banked undercooling without discarding the cells or the solute
field that took ~1e6 steps to grow.

FILE FORMAT (verified against the solver source and by byte count against two
real frames).  Despite the .vtk extension these files are not VTK:

    line 1 : nx            ASCII.  The solver SKIPS all three lines
    line 2 : ny                    (fscanf "%*[^\\n]\\n" x3) and never
    line 3 : nz                    validates them against the mesh.
    phi[0] : nx*ny*nz float64 big-endian, then a single b'\\n'
    phi[1] : same                                       + b'\\n'
    mu[0]  : same   (NUMCOMPONENTS-1 blocks)            + b'\\n'
    c[0]   : same   (NUMCOMPONENTS-1 blocks)            + b'\\n'
    T      : same   (only when !ISOTHERMAL), NO trailing separator

Element order is x outer, z middle, y inner:  index = y + ny*(z + nz*x).
Reader: solverloop/file_writer.h:521 read_cells_vtk_mpi_binary.

Byte-count check, 2 phases / 2 components / non-isothermal -> 5 blocks:
    5*8*nx*ny*nz + len(header) + 4 separators
    3168x2240x1 -> 283,852,816 (matches stage-5 frame exactly)
     768x5600x1 -> 172,032,015 (matches pilot frame 0 exactly)

The reader has no header check, but the b'\\n' separators are a real integrity
check: a wrong mesh gives "Missing binary separator", a short file gives
"Failed to read binary".  Mistakes fail loudly.
"""

import argparse
import os
import sys

import numpy as np

BE = np.dtype(">f8")
SEP = b"\n"


# --------------------------------------------------------------------------
# I/O
# --------------------------------------------------------------------------
def read_frame(path, nphases=2, ncomp_minus1=1, has_temperature=True):
    """Read a MicroSim restart frame.  Returns (meta, fields)."""
    size = os.path.getsize(path)
    with open(path, "rb") as fp:
        header = [fp.readline() for _ in range(3)]
        nx, ny, nz = (int(line) for line in header)
        header_bytes = sum(len(line) for line in header)

        nblocks = nphases + 2 * ncomp_minus1 + (1 if has_temperature else 0)
        nsep = nphases + 2 * ncomp_minus1  # temperature has no separator
        n = nx * ny * nz
        expect = 8 * n * nblocks + header_bytes + nsep
        if size != expect:
            raise ValueError(
                f"{path}: size {size} != expected {expect} for "
                f"{nx}x{ny}x{nz}, {nblocks} blocks.  Wrong NUMPHASES / "
                f"NUMCOMPONENTS / ISOTHERMAL?"
            )

        def block(name, sep=True):
            arr = np.fromfile(fp, dtype=BE, count=n)
            if arr.size != n:
                raise ValueError(f"{path}: short read in {name}")
            if sep:
                got = fp.read(1)
                if got != SEP:
                    raise ValueError(
                        f"{path}: missing separator after {name}, got {got!r}"
                    )
            # (nx, ny, nz) with y fastest -> reshape (nx, nz, ny) then swap
            return arr.astype(np.float64).reshape(nx, nz, ny).transpose(0, 2, 1)

        fields = {}
        for a in range(nphases):
            fields[f"phi{a}"] = block(f"phi[{a}]")
        for k in range(ncomp_minus1):
            fields[f"mu{k}"] = block(f"mu[{k}]")
        for k in range(ncomp_minus1):
            fields[f"c{k}"] = block(f"c[{k}]")
        if has_temperature:
            fields["T"] = block("T", sep=False)

    meta = dict(nx=nx, ny=ny, nz=nz, nphases=nphases,
                ncomp_minus1=ncomp_minus1, has_temperature=has_temperature,
                path=path, size=size)
    return meta, fields


def write_frame(path, fields, meta):
    """Write a MicroSim restart frame in the byte layout read_frame expects."""
    nx, ny, nz = meta["nx"], meta["ny"], meta["nz"]
    order = ([f"phi{a}" for a in range(meta["nphases"])]
             + [f"mu{k}" for k in range(meta["ncomp_minus1"])]
             + [f"c{k}" for k in range(meta["ncomp_minus1"])])
    tmp = path + ".partial"
    with open(tmp, "wb") as fp:
        fp.write(f"{nx}\n{ny}\n{nz}\n".encode())
        for name in order:
            arr = fields[name]
            if arr.shape != (nx, ny, nz):
                raise ValueError(f"{name} shape {arr.shape} != {(nx, ny, nz)}")
            arr.transpose(0, 2, 1).astype(BE).tofile(fp)
            fp.write(SEP)
        if meta["has_temperature"]:
            arr = fields["T"]
            arr.transpose(0, 2, 1).astype(BE).tofile(fp)
    os.replace(tmp, path)
    return os.path.getsize(path)


# --------------------------------------------------------------------------
# geometry
# --------------------------------------------------------------------------
def interface_row(phi_solid, threshold=0.5):
    """Topmost row per column where the solid phase field exceeds threshold."""
    solid = phi_solid[:, :, 0] > threshold           # (nx, ny)
    any_solid = solid.any(axis=1)
    top = np.where(any_solid, solid.shape[1] - 1 - solid[:, ::-1].argmax(axis=1), 0)
    return top.astype(np.float64)


def best_crop_offset(fields, width, report=8):
    """
    Choose the x offset whose wrap seam is least visible.

    Cropping columns [x0, x0+W) and re-imposing periodicity makes column
    x0+W-1 adjacent to x0.  In the source its neighbour was column x0+W.  So
    the seam mismatch is the distance between column x0 and column x0+W --
    that is the column the wrap replaces.  Scored against the mean distance
    between genuinely adjacent columns, so a ratio near 1.0 means the seam is
    indistinguishable from an ordinary neighbour step.
    """
    phi = fields["phi0"][:, :, 0]
    c = fields["c0"][:, :, 0]
    nx = phi.shape[0]
    c_scale = c.std() or 1.0

    def dist(i, j):
        return (np.sqrt(np.mean((phi[i] - phi[j]) ** 2))
                + np.sqrt(np.mean((c[i] - c[j]) ** 2)) / c_scale)

    idx = np.arange(nx)
    baseline = np.mean([dist(i, (i + 1) % nx) for i in idx])
    cost = np.array([dist(x0, (x0 + width) % nx) for x0 in idx])

    best = int(np.argmin(cost))
    ranking = np.argsort(cost)[:report]
    return best, cost, baseline, ranking


# --------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("source", help="restart frame to re-embed")
    p.add_argument("--out", required=True, help="output frame (name it <prefix>_0.vtk)")
    p.add_argument("--nx", type=int, required=True, help="cropped width in cells")
    p.add_argument("--ny", type=int, required=True, help="new domain height in cells")
    p.add_argument("--x-offset", type=int, default=None,
                   help="force a crop offset instead of the seam-matched one")
    p.add_argument("--pad-mode", choices=("replicate", "c0"), default="replicate",
                   help="replicate: copy the source's top liquid row (mu and c "
                        "mutually consistent by construction).  c0: pad at the "
                        "nominal alloy composition, with mu interpolated from "
                        "the source's own liquid cells.")
    p.add_argument("--c0", type=float, default=0.052953,
                   help="nominal Si mole fraction, for --pad-mode c0")
    p.add_argument("--base-temp", type=float, default=884.367547,
                   help="only used to write a self-consistent T block; the "
                        "solver overwrites it from Tempgrady")
    p.add_argument("--gradient", type=float, default=1.1e6, help="K/m")
    p.add_argument("--dx", type=float, default=1.0771e-8, help="m")
    p.add_argument("--self-test", action="store_true",
                   help="round-trip the source with no crop and no pad and "
                        "byte-compare against the original")
    args = p.parse_args()

    meta, fields = read_frame(args.source)
    nx0, ny0, nz = meta["nx"], meta["ny"], meta["nz"]
    print(f"source     {args.source}")
    print(f"  mesh     {nx0} x {ny0} x {nz}   ({meta['size']:,} bytes, layout OK)")

    if args.self_test:
        out = args.out
        n = write_frame(out, fields, meta)
        same = (open(args.source, "rb").read() == open(out, "rb").read())
        print(f"SELF-TEST  rewrote {n:,} bytes; byte-identical to source: {same}")
        sys.exit(0 if same else 1)

    if args.ny < ny0:
        sys.exit(f"ERROR: --ny {args.ny} < source ny {ny0}; this tool pads, "
                 f"it does not truncate.")
    if args.nx > nx0:
        sys.exit(f"ERROR: --nx {args.nx} > source nx {nx0}")

    top = interface_row(fields["phi0"])
    print(f"  tips     row {top.max():.0f}   roots row {top.min():.0f}   "
          f"mean {top.mean():.1f}  (groove depth {top.max()-top.min():.0f} cells)")

    # ---- crop -------------------------------------------------------------
    if args.nx == nx0:
        x0 = 0
        print("  crop     none (full width kept)")
    else:
        if args.x_offset is None:
            x0, cost, baseline, ranking = best_crop_offset(fields, args.nx)
            print(f"  seam     best x-offset {x0}: cost {cost[x0]:.4g} vs mean "
                  f"adjacent-column cost {baseline:.4g}  "
                  f"(ratio {cost[x0]/baseline:.2f})")
            print(f"           worst offset would be {cost.max()/baseline:.1f}x; "
                  f"next best offsets {list(ranking[1:5])}")
        else:
            x0 = args.x_offset
            _, cost, baseline, _ = best_crop_offset(fields, args.nx)
            print(f"  seam     forced x-offset {x0}: ratio "
                  f"{cost[x0]/baseline:.2f}")
        idx = (np.arange(args.nx) + x0) % nx0
        fields = {k: v[idx] for k, v in fields.items()}

    # ---- pad --------------------------------------------------------------
    npad = args.ny - ny0
    out = {}
    if args.pad_mode == "replicate":
        pad_from = ny0 - 1
        pad_vals = {k: fields[k][:, pad_from:pad_from + 1, :] for k in fields}
        src_c = float(np.mean(fields["c0"][:, pad_from, 0]))
        print(f"  pad      {npad} rows by replicating source row {pad_from} "
              f"(mean x_Si {src_c:.6f})")
    else:
        liquid = fields["phi1"][:, :, 0] > 0.99
        cl = fields["c0"][:, :, 0][liquid]
        mul = fields["mu0"][:, :, 0][liquid]
        srt = np.argsort(cl)
        mu_at_c0 = float(np.interp(args.c0, cl[srt], mul[srt]))
        near = np.abs(cl - args.c0).min()
        print(f"  pad      {npad} rows at x_Si = {args.c0:.6f}, "
              f"mu interpolated from {liquid.sum():,} source liquid cells "
              f"-> {mu_at_c0:.8e}  (nearest source cell differs by {near:.2e})")
        pad_vals = {}
        for k in fields:
            v = np.zeros((args.nx, 1, nz))
            if k == "phi1":
                v[:] = 1.0
            elif k == "c0":
                v[:] = args.c0
            elif k == "mu0":
                v[:] = mu_at_c0
            elif k == "T":
                v[:] = fields["T"][:, -1:, :]
            pad_vals[k] = v

    for k, v in fields.items():
        out[k] = np.concatenate([v, np.repeat(pad_vals[k], npad, axis=1)], axis=1)

    # ---- rewrite T self-consistently (the solver overwrites it anyway) -----
    y = np.arange(args.ny)
    out["T"] = np.broadcast_to(
        (args.base_temp + args.gradient * y * args.dx)[None, :, None],
        (args.nx, args.ny, nz)).copy()

    meta_out = dict(meta, nx=args.nx, ny=args.ny)
    n = write_frame(args.out, out, meta_out)
    print(f"wrote      {args.out}")
    print(f"  mesh     {args.nx} x {args.ny} x {nz}   ({n:,} bytes)")

    # ---- verify by reading our own output back ----------------------------
    m2, f2 = read_frame(args.out)
    top2 = interface_row(f2["phi0"])
    phisum = float(np.abs(f2["phi0"] + f2["phi1"] - 1.0).max())
    print(f"  reread   OK; tips row {top2.max():.0f}; "
          f"max |sum(phi)-1| = {phisum:.3e}")
    print(f"  headroom {args.ny - top2.max():.0f} cells above the tips = "
          f"{(args.ny - top2.max()) * args.gradient * args.dx:.2f} K")


if __name__ == "__main__":
    main()
