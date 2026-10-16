"""红旗核查：impact_speed / momentum 是否随 h 按 sqrt(2gh) 走？"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

from climbing.coupling.opensim_grf import ground_reaction  # noqa: E402

print("  h(m)   v_impact   sqrt(2gh)   ratio   momentum(N·s)  m*sqrt(2gh)  ratio   GRFtot(kN)")
for h in (0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 8.0, 10.0):
    g = ground_reaction(height_m=h)
    v = float(g.impact_speed_ms)
    v0 = float(np.sqrt(2 * 9.81 * h))
    m0 = g.mass_kg * v0
    print(f"  {h:4.1f}  {v:9.4f}  {v0:10.4f}  {v/v0:6.3f}  "
          f"{g.momentum_ns:12.2f}  {m0:11.2f}  {g.momentum_ns/m0:6.3f}  "
          f"{g.peak_total_n/1e3:9.2f}")
print()
g = ground_reaction(height_m=2.0)
print("h=2 全字段：")
for f in ("t_s", "f_total_n", "f_left_n", "f_right_n", "height_m", "mass_kg",
          "posture", "on_pad", "impact_speed_ms", "peak_total_n", "impulse_ns",
          "momentum_ns", "impulse_rel_err", "window_s", "split"):
    v = getattr(g, f)
    if isinstance(v, np.ndarray):
        print(f"  {f}: ndarray shape={v.shape} max={v.max():.4g} t0={g.t_s[0]} t1={g.t_s[-1]}")
    else:
        print(f"  {f}: {v}")
