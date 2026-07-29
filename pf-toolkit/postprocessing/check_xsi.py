import sys, numpy as np
sys.path.insert(0, sys.argv[2])
from pathlib import Path
from analyze_pilot_morphology import read_fragment
dims, arr = read_fragment(Path(sys.argv[1]))
x = arr["Composition_Si"][:, :, 0]
phi = arr["FCC_A1"][:, :, 0]
n = x.size
X_EUT = 0.12074      # mol frac Si at eutectic (12.506 wt%)
for lab, m in [("xSi < 0", x < 0), ("xSi < -1e-3", x < -1e-3),
               ("xSi > x_eut", x > X_EUT)]:
    c = int(m.sum())
    print(f"{lab:15s} {c:9d} cells  {100*c/n:7.4f}%")
print(f"min {x.min():.6f}  max {x.max():.6f}  mean {x.mean():.6f}")
neg = x < 0
if neg.any():
    print(f"phi at negative cells: min {phi[neg].min():.4f} max {phi[neg].max():.4f} mean {phi[neg].mean():.4f}")
    ys = np.where(neg)[1]
    print(f"y-range of negative cells: {ys.min()} .. {ys.max()} cells")
hi = x > X_EUT
if hi.any():
    ys = np.where(hi)[1]
    print(f"y-range of >x_eut cells: {ys.min()} .. {ys.max()} cells; phi mean {phi[hi].mean():.4f}")
