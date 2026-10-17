"""★ CANONICAL 复算：首折 `h*`（**契约正确**的表格式）。

契约（`CANONICAL.md` §0）：`m_theta_table = [[θ_deg, M_Nm], ...]`
—— 2 维、两列、θ 用**度数**、严格升序。（`pad.py` 自己 `np.radians(_mt[:, 0])`。）

自检（两道，缺一不可）：
  ① 正面：复刻 `pad.py` 的转换后 `interp(1.0 rad)` 必须 ≈ 254.7 N·m；
  ② 反面：传 `(x_rad, y)` 元组必须**抛 ValueError**（硬校验生效）。

产物：`results/hstar_canonical.json`
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

ROOT = Path(r'D:\Project\climbing_fall_analysis')
sys.path.insert(0, str(ROOT / "src"))
import climbing.pad as P  # noqa: E402

MASS, L, V_REF, M_ECC, K_TH, THR = 80.0, 0.425, 1.0, 540.0, 350.0, 6200.0
ANCH = [(0.0, 1.30), (15.0, 1.79), (45.0, 2.90), (60.0, 3.20), (90.0, 2.78), (135.0, 2.78)]
_DEG = np.arange(0.0, 135.0 + 1e-9, 2.5)
_FIT = PchipInterpolator([a for a, _ in ANCH], [b * MASS for _, b in ANCH])
#: ★ 契约格式：`[[θ_deg, M_Nm], ...]`
TABLE = [[float(t), float(_FIT(t))] for t in _DEG]
COM = {"soft": 0.43, "self": 0.37, "stiff": 0.23}
KEYS = ("k_theta_nm_per_rad", "segment_l_m", "theta_max_rad",
        "m_ecc_max_nm", "ecc_rate_ref_ms", "m_theta_table")
GRID = (0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0)

# ---- ① 正面自检：复刻 pad.py 的转换 -------------------------------------
_arr = np.asarray(TABLE, dtype=float)
_v100 = float(np.interp(1.0, np.radians(_arr[:, 0]), _arr[:, 1]))
assert abs(_v100 - 254.7) < 3.0, f"① 表自检失败：interp(1.0 rad) = {_v100}（期望 ~254.7）"
print(f"① 正面自检：interp(1.0 rad) = {_v100:.1f} N·m（期望 ~254.7）OK")


def run(lab, h, table=None):
    prev = {k: P.TWO_SEG_LEG.get(k) for k in KEYS}
    c0 = P.RIGID_CONTACT.get("k_n_per_m")
    try:
        P.TWO_SEG_LEG.update({
            "k_theta_nm_per_rad": K_TH, "segment_l_m": L,
            "theta_max_rad": float(2.0 * np.arccos(1.0 - COM[lab] / (2.0 * L))),
            "m_ecc_max_nm": M_ECC, "ecc_rate_ref_ms": V_REF,
            "m_theta_table": TABLE if table is None else table})
        P.RIGID_CONTACT["k_n_per_m"] = 2.1e6
        return P.simulate_boulder_fall(height_m=h, mass_kg=MASS, on_pad=False,
                                       posture=f"anchored-{lab}")
    finally:
        for k in KEYS:
            P.TWO_SEG_LEG[k] = prev[k]
        P.RIGID_CONTACT["k_n_per_m"] = c0


# ---- ② 反面自检：错格式必须被拒 ----------------------------------------
try:
    run("self", 0.5, table=(np.radians(_arr[:, 0]), _arr[:, 1]))
    print("② 反面自检：❌ 未抛错 —— 硬校验失效（传元组会被静默读成废表）！")
    raise SystemExit(1)
except ValueError as e:
    print(f"② 反面自检：OK —— 错格式被拒（{str(e)[:52]}…）")

# ---- 主算：细网格 -------------------------------------------------------
print()
print(f"  {'lab':>6} " + "".join(f"{h:>8.2f}" for h in GRID) + f"{'h* (m)':>12}")
out = {"threshold_n": THR, "grid_m": list(GRID), "com_depth_m": COM,
       "table_contract": "[[θ_deg, M_Nm], ...]", "table_interp_1rad_nm": _v100}
for lab in COM:
    row = {h: float(run(lab, h).meta["f_flex_peak_n"]) for h in GRID}
    hs = sorted(row)
    hsx = None
    for i in range(len(hs) - 1):
        a, b = row[hs[i]], row[hs[i + 1]]
        if a <= THR <= b:
            hsx = hs[i] + (THR - a) * (hs[i + 1] - hs[i]) / (b - a)
            break
    out[lab] = {str(h): row[h] for h in hs}
    out[lab + "__hstar_m"] = hsx
    print(f"  {lab:>6} " + "".join(f"{row[h] / 1e3:>8.3f}" for h in hs)
          + f"{('%.3f' % hsx) if hsx else '< 0.05':>12}")
print("\n  （单位 kN；判据 = 跟骨阈 6200 N）")
dst = ROOT / "results" / "hstar_canonical.json"
dst.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
print(f"\nOUT = {dst.relative_to(ROOT)}")
