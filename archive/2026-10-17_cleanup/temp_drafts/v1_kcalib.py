"""★★ 用文献标定 `k_flex`（用户要求「有没有文献可以标定真实刚度」）。

标定锚（`PMC9465758`，You et al.，24 名男性长曲棍球运动员，**0.42 m 落箱**，双臂落地）：
  | 落地风格 | 腿刚度 (kN/m) | 触地时间 (s) | COM 深度 (m) | **GRFmax (N/kg)** |
  | Soft         |  82.4 ± 25.3 | 0.52 | 0.43 | **33.13** (≈3.38 BW) |
  | Self-selected| 107.7 ± 39.9 | 0.45 | 0.37 | 36.26 (≈3.70 BW) |
  | **Stiff**    | **243.3 ± 55.3** | **0.29** | **0.23** | **52.55** (≈**5.36 BW**) |

另：DeVita & Skelly 1992（59 cm 落，soft vs stiff）——
  **stiff 落地 GRFmax ≈ 32 N/kg ≈ 3.26 BW，soft ≈ 23 N/kg ≈ 2.35 BW**（50 ms 处）。
另：Arampatzis 等（跑步）leg stiffness **25.3–35.2 kN/m**、vertical stiffness **30.8–93.0 kN/m**；
  短跑 vertical stiffness **74–171 kN/m**。

本件：把我们的 `pad.py` 放在 **同高度（0.42 m）**，比 GRFmax/BW，并给出「模型 k vs 文献 k」对照。
"""
import sys
from pathlib import Path
from dataclasses import replace

ROOT = Path(r'D:\Project\climbing_fall_analysis')
sys.path.insert(0, str(ROOT / "src"))
import climbing.pad as P  # noqa: E402

M, G = 80.0, 9.80665
BW = M * G

print('=== ① 模型 `k_flex` vs 文献腿刚度 ===\n')
LIT = [("Soft landing", 82.4e3), ("Self-selected", 107.7e3), ("★ Stiff landing", 243.3e3),
       ("Stiff +1SD", 298.6e3), ("跑步 leg stiffness 25–35 kN/m", 30.0e3),
       ("短跑 vertical stiffness 74–171 kN/m", 120.0e3)]
MODEL = [("controlled-drop", P.POSTURES["controlled-drop"].k_flex),
         ("feet-first-stiff", P.POSTURES["feet-first-stiff"].k_flex)]
print(f'  {"来源":<34} {"k (N/m)":>12} {"k (kN/m)":>10}')
for nm, k in LIT:
    print(f'  {nm:<34} {k:>12.1e} {k/1e3:>10.1f}')
print('  ' + '-' * 58)
for nm, k in MODEL:
    print(f'  【模型】{nm:<27} {k:>12.1e} {k/1e3:>10.1f}')
print()
for nm, k in MODEL:
    r = k / 243.3e3
    print(f'  ⇒ 模型 {nm} 是实测**刚性落地**（243.3 kN/m）的 **{r:.2f}×**'
          f'；是短跑 vertical 上限（171）的 {k/171e3:.2f}×')

print('\n=== ② 同高度对照（0.42 m = 文献落箱高度）===\n')
print(f'  文献：soft GRFmax=3.38 BW；stiff GRFmax=5.36 BW（DeVita59cm: 3.26/2.35 BW）')
print(f'\n  {"姿势":<22} {"k(kN/m)":>9} {"F_pad(kN)":>10} {"F_pad/BW":>10} '
      f'{"F_flex/BW":>11} {"vs 实测 stiff":>13}')
for nm in ("controlled-drop", "feet-first-stiff"):
    r = P.simulate_boulder_fall(height_m=0.42, mass_kg=M, on_pad=False, posture=nm)
    fp = float(r.peak_force_n)
    ff = float(r.meta.get("f_flex_peak_n", 0.0))
    k = P.POSTURES[nm].k_flex
    print(f'  {nm:<22} {k/1e3:>9.0f} {fp/1e3:>10.2f} {fp/BW:>10.2f} {ff/BW:>11.2f} '
          f'{(fp/BW)/5.36:>12.2f}×')

print('\n=== ③ 若把 `k_flex` 设到文献值（stiff 243.3 kN/m），我们的 1 m 坠落给多少？===\n')
base = P.POSTURES["feet-first-stiff"]
for nm, k in (("文献 soft 82.4", 82.4e3), ("文献 self 107.7", 107.7e3),
              ("文献 stiff 243.3", 243.3e3)):
    key = f"lit_{nm}"
    P.POSTURES[key] = replace(base, name=nm, k_flex=k, x_flex_max=0.6)
    r = P.simulate_boulder_fall(height_m=1.0, mass_kg=M, on_pad=False, posture=key)
    fp, ff = float(r.peak_force_n), float(r.meta.get("f_flex_peak_n", 0.0))
    print(f'  k={nm:<18} F_pad={fp/1e3:>7.2f} kN ({fp/BW:>5.1f} BW)   '
          f'F_flex={ff/1e3:>7.2f} kN   ⇒ 跟骨阈 6.2 kN 的 {ff/6200:.2f}×')
print('\n参照：目标带 6–9 m；尸体跟骨阈 6.2 kN')
