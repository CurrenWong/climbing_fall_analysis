"""★ V1 的 **1D 力判据近似**：首折 = h 使「单侧跟骨力」达到文献阈值。

语义（`src/climbing/bone.py`）：
  * `calcaneus.load_fraction = 0.50`（双足着地，单侧跟骨承担约一半）
  * 跟骨力 = `0.50 × F_pad(h)`（`f_pad` = 整通路总力）
  * **名义阈值** = σ_c × A = 150 MPa × 900 mm² = **135 kN**（`knockdown_report` 自评高估 **16.9–33.8×**）
  * **文献阈值**（Barnes IRCOBI 2019 / Voo HIPC）= **4–8 kN** ⇒ **可辩护阈值**

本件算两条：
  (a) 名义阈值路线（项目 `assess_sites` 现用）→ 首折
  (b) **文献阈值路线（可辩护的 1D 近似）** → 首折
"""
import inspect
import sys
from pathlib import Path

ROOT = Path(r'D:\Project\climbing_fall_analysis')
sys.path.insert(0, str(ROOT / "src"))

import climbing.pad as P  # noqa: E402
import climbing.bone as B  # noqa: E402

print('simulate_boulder_fall 签名：', inspect.signature(P.simulate_boulder_fall))
site = B.BONE_SITES["calcaneus"]
print(f'\n跟骨 site：A={site.area_mm2} mm²  σ_c={site.sigma_mpa} MPa  '
      f'load_fraction={site.load_fraction}  confidence={site.confidence}')
print(f'  名义阈值 = {site.force_threshold_n/1e3:.1f} kN')
print(f'  ⇒ 需要的总 GRF（名义）= {site.force_threshold_n/site.load_fraction/1e3:.1f} kN')
print(f'  文献阈值 4–8 kN ⇒ 需要的总 GRF = {4e3/site.load_fraction/1e3:.1f}–'
      f'{8e3/site.load_fraction/1e3:.1f} kN')

# --- 取 F_pad(h) ---
hs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15]
rows = []
for h in hs:
    try:
        r = P.simulate_boulder_fall(height_m=float(h))
    except TypeError:
        r = P.simulate_boulder_fall(float(h))
    # 找 peak_pad_n
    pk = None
    if hasattr(r, "peak_pad_n"):
        pk = r.peak_pad_n
    elif isinstance(r, dict):
        pk = r.get("peak_pad_n") or (r.get("meta") or {}).get("peak_pad_n")
    else:
        meta = getattr(r, "meta", None)
        if isinstance(meta, dict):
            pk = meta.get("peak_pad_n")
    if pk is None:
        print('  [warn] 取不到 peak_pad_n，r 类型 =', type(r))
        print('  dir:', [a for a in dir(r) if not a.startswith('_')][:30])
        break
    rows.append({"h": h, "f_pad_n": float(pk), "calc_n": float(pk) * site.load_fraction})

print(f'\n  {"h(m)":>5} {"F_pad(kN)":>11} {"单侧跟骨(kN)":>13} {"名义判定":>9} {"文献判定":>9}')
for r in rows:
    cn = r["calc_n"]
    print(f'  {r["h"]:>5} {r["f_pad_n"]/1e3:>11.2f} {cn/1e3:>13.2f} '
          f'{"折" if cn >= site.force_threshold_n else "—":>9} '
          f'{"折(≤4kN)" if cn >= 4e3 else "—":>9}')


def cross(thr_n):
    o = sorted(rows, key=lambda r: r["h"])
    for a, b in zip(o, o[1:]):
        if a["calc_n"] < thr_n <= b["calc_n"]:
            return a["h"] + (b["h"] - a["h"]) * (thr_n - a["calc_n"]) / (b["calc_n"] - a["calc_n"])
    return None


print('\n=== 首折（1D 力判据）===')
for nm, thr in (("名义 135 kN（项目 assess_sites 现用）", site.force_threshold_n),
                ("文献上界 8 kN", 8e3), ("文献下界 4 kN", 4e3)):
    h = cross(thr)
    print(f'  {nm:<36} ⇒ 首折 = ' + ('未跨' if h is None else f'{h:.2f} m'))
print('\n参照：FE 单骨线（最可辩护上界）≲3.05 m；目标带 6–9 m')
