"""P1 精测：F_subt 与 GRF、每脚 GRF、体重的关系。"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

from climbing.coupling.joint_loads import subtalar_reaction      # noqa: E402
from climbing.coupling.opensim_fall import run_dead_drop          # noqa: E402
from climbing.coupling.opensim_grf import ground_reaction         # noqa: E402

print("  h(m)  GRFtot(kN)  GRF_r(kN)  F_subt(kN)  F/GRFtot  F/GRF_r  F/BW  split")
rows = []
for h in (1.0, 2.0, 3.0, 5.0, 10.0):
    grf = ground_reaction(height_m=h)
    fall = run_dead_drop(grf)
    jr = subtalar_reaction(fall, grf, side="r")
    f = float(jr.peak_vertical_n)
    gt = float(grf.peak_total_n)
    gr = float(np.max(grf.f_right_n))
    bw = grf.mass_kg * 9.81
    rows.append((h, gt, gr, f, f / gt, f / gr, f / bw))
    print(f"  {h:4.1f}  {gt/1e3:10.2f}  {gr/1e3:9.2f}  {f/1e3:10.2f}  "
          f"{f/gt:8.3f}  {f/gr:7.3f}  {f/bw:4.1f}  {grf.split}"
          f"  loaded={grf.loaded_sides}  share={grf.peak_share}")

print()
print("质量 =", grf.mass_kg, "kg  ⇒ BW =", grf.mass_kg * 9.81, "N")
print("冲击速度 h=2:", grf.impact_speed_ms, "m/s")
print("momentum =", grf.momentum_ns, "N·s")
print()
print("判读：F/GRFtot ≈ 1.0 ⇒ 单脚承全量（右脚独担）；≈0.5 ⇒ 双足分担，每脚一半")
