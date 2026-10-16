"""B1 预演：把 `x_flex_max` 锚定到**实测 COM 深度**（PMC9465758：soft 0.43 / self 0.37 / stiff 0.23 m）
后，算力判据首折 h*；并**标出止段（`xf > x_flex_max`）被触发的格子**（那些格不可用）。
只读 + 临时姿势对象/临时门（finally 恢复）。**不改仓库。**"""
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

ROOT = Path(r'D:\Project\climbing_fall_analysis')
sys.path.insert(0, str(ROOT / "src"))
import climbing.pad as P  # noqa: E402

G, LEG0 = 9.80665, 0.85
BW = 80.0 * G
F_THR = 6200.0
GRID = [0.20, 0.50, 1.00, 2.00, 5.00, 7.00, 9.00]
# 实测 COM 深度（`PMC9465758`「Drop jump depth」列）
COM = {"soft": 0.43, "self": 0.37, "stiff": 0.23}
LAB = ("soft", "self", "stiff")


def run(post, h, kn):
    o = P.RIGID_CONTACT["k_n_per_m"]
    try:
        P.RIGID_CONTACT["k_n_per_m"] = kn
        return P.simulate_boulder_fall(height_m=h, mass_kg=80.0, on_pad=False, posture=post)
    finally:
        P.RIGID_CONTACT["k_n_per_m"] = o


def stop_engaged(r, xm):
    """冲击窗内 `xf` 是否超过 `x_flex_max`（= 进 k_stop 止段）。"""
    t = np.asarray(r.t_s, float)
    zf, zt = np.asarray(r.z_foot_m, float), np.asarray(r.z_torso_m, float)
    t0, t1 = r.meta["impact_window_s"]
    i0, i1 = int(np.searchsorted(t, t0)), max(int(np.searchsorted(t, t1)), 2)
    xf = np.maximum(0.0, LEG0 - (zf - zt))
    return float(xf[i0:i1].max()), bool(xf[i0:i1].max() > xm)


def hstar(vals):
    for i in range(len(GRID) - 1):
        if vals[i] <= F_THR <= vals[i + 1]:
            return GRID[i] + (F_THR - vals[i]) / (vals[i + 1] - vals[i]) * (GRID[i + 1] - GRID[i])
    return None


posts = {l: replace(P.POSTURES[f"anchored-{l}"], name=f"B1-{l}",
                    x_flex_max=COM[l]) for l in LAB}
print("=== B1：x_flex_max 锚定到实测 COM 深度 ===")
for l in LAB:
    print(f"  {l:<6} k_flex={posts[l].k_flex/1e3:>5.1f} kN/m   "
          f"x_flex_max={posts[l].x_flex_max:.2f} m（实测 COM 深度）")
print(f"  （k_stop 仍为 pad.py 硬编码 3.0e5 N/m ⇒ 标出被触发的格子）")

for kn in (None, 2.1e6):
    tag = "刚性接触（= LEG_K_ANCHOR 口径）" if kn is None else f"实测跟垫接触 k_n={kn:.1e}"
    print(f"\n=== {tag} ===")
    print(f"  {'档':<6} " + " ".join(f"{h:>7.2f}m" for h in GRID) + f" {'h*(m)':>8} {'止段?':>16}")
    for l in LAB:
        vals, eng = [], []
        for h in GRID:
            r = run(posts[l], h, kn)
            vals.append(float(r.meta["f_flex_peak_n"]))
            xfm, e = stop_engaged(r, posts[l].x_flex_max)
            eng.append(e)
        hs = hstar(vals)
        # 单调性 + 首个止段高度
        mono = all(vals[i] < vals[i + 1] for i in range(len(vals) - 1))
        first_eng = next((f"{GRID[i]:.2f}" for i, e in enumerate(eng) if e), "无")
        print(f"  {l:<6} " + " ".join(f"{v/1e3:>7.2f} " for v in vals)
              + f" {('%.3f' % hs) if hs else '>9':>8} "
              + f" {'首止段@' + first_eng + 'm':>16}"
              + ("" if mono else "  ⚠非单调"))
    print(f"  （止段格 = `xf > x_flex_max` ⇒ 峰值被 k_stop 顶起，**该格不可用**）")

print("\n=== 参照 ===")
print("  阈值 6.2 kN（尸体跟骨，Yoganandan）｜流行病学 1.1–4 m｜临床 ≤5 ft(1.52 m) 亦可致骨折｜[Y25] 靶 7~9 m")
