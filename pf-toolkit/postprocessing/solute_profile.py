#!/usr/bin/env python3
"""Horizontally-averaged Si profile vs y, to measure whether the moving-window
refill is truncating the solute boundary layer.

The mass audit in CL_Shift.h cannot detect this: expected_after is DEFINED to
include refill at c0, so its residual is ~1e-13 whether or not real enrichment
is being deleted. This measures the physical quantity instead.
"""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(sys.argv[0]).resolve().parent))
from analyze_pilot_morphology import read_fragment

C0 = 0.052953          # cfill liquid, mol frac Si
DX = 1.0771e-8
LD_CELLS = 6.25e-9 / 1.3333e-2 / DX      # D/V in cells

frame = Path(sys.argv[1])
shiftj = int(sys.argv[2]) if len(sys.argv) > 2 else 1648

dims, arr = read_fragment(frame)
# Array shape is (nx, ny, nz) with nz = 1 for these 2-D runs; axis 1 is y.
if dims[2] != 1:
    raise SystemExit("2-D frames only")
c = arr["Composition_Si"][:, :, 0]
sol = arr["FCC_A1"][:, :, 0]
cy = c.mean(axis=0)
soly = sol.mean(axis=0)
ny = cy.size

# tip = highest row containing any solid
solid_rows = np.where(sol.max(axis=0) > 0.5)[0]
tip = int(solid_rows.max()) if solid_rows.size else -1

print(f"{frame.name}  ny={ny}  tip_row={tip}  Shiftj={shiftj}")
print(f"  l_D = {LD_CELLS:.1f} cells; headroom above Shiftj = {ny-shiftj} cells "
      f"= {(ny-shiftj)/LD_CELLS:.2f} l_D")
print(f"  liquid-only mean c in top 20 rows = {cy[-20:].mean():.8f} "
      f"({100*(cy[-20:].mean()/C0-1):+.4f}% vs c0)")
print()
print("   row   dist_above_tip[cells]  [l_D]   <c>/c0-1 [%]   <phi_FCC>")
for row in [ny-1, ny-5, ny-20, ny-60, ny-120, ny-180, ny-240, ny-300, ny-400]:
    if row <= 0 or row >= ny:
        continue
    d = row - tip
    print(f"  {row:5d} {d:18d} {d/LD_CELLS:8.2f} {100*(cy[row]/C0-1):+13.4f}"
          f" {soly[row]:12.3e}")

# Where has the layer decayed to within 0.1% of c0?
above = np.arange(ny) > tip
excess = np.abs(cy/C0 - 1)
ok = np.where(above & (excess < 1e-3))[0]
if ok.size:
    print(f"\n  first row above tip within 0.1% of c0: {ok.min()} "
          f"= {(ok.min()-tip)/LD_CELLS:.2f} l_D above the tip")
else:
    print("\n  NO row above the tip is within 0.1% of c0")
print(f"  excess at the very top row  : {100*(cy[-1]/C0-1):+.4f}%")
print(f"  excess at Shiftj row ({shiftj}): {100*(cy[shiftj]/C0-1):+.4f}%"
      if shiftj < ny else "")
