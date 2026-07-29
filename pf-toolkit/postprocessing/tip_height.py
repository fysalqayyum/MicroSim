import sys, numpy as np
from pathlib import Path
sys.path.insert(0, sys.argv[1])
from analyze_pilot_morphology import read_fragment, timestep
DX_UM, DT = 0.010771, 2.0e-9
prev = None
print(f"{'step':>9} {'h_max':>8} {'h_mean':>8} {'h_min':>8} {'v_tip':>9} {'v_mean':>9}")
print(f"{'':9} {'cells':>8} {'cells':>8} {'cells':>8} {'mm/s':>9} {'mm/s':>9}")
for p in sorted((Path(sys.argv[2]).glob("*.vtk")), key=timestep):
    st = timestep(p)
    if st < int(sys.argv[3]):
        continue
    dims, arr = read_fragment(p)
    phi = arr["FCC_A1"][:, :, 0]           # [x, y]
    solid = phi > 0.5
    # highest solid row per column; columns with no solid -> 0
    hmax_col = np.where(solid.any(axis=1), solid.shape[1] - 1 - np.argmax(solid[:, ::-1], axis=1), 0)
    hmax, hmean, hmin = hmax_col.max(), hmax_col.mean(), hmax_col.min()
    if prev is None:
        vt = vm = float("nan")
    else:
        dt_s = (st - prev[0]) * DT
        vt = (hmax - prev[1]) * DX_UM * 1e-6 / dt_s * 1e3
        vm = (hmean - prev[2]) * DX_UM * 1e-6 / dt_s * 1e3
    print(f"{st:9d} {hmax:8.0f} {hmean:8.1f} {hmin:8.0f} {vt:9.3f} {vm:9.3f}")
    prev = (st, hmax, hmean)
