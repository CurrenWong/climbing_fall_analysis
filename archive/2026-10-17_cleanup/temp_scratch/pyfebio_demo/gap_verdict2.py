"""独立裁决 v2：顶面相对距骨的 z 位移方向 —— 不依赖任何符号约定。"""
import os
import re
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
PY = ROOT / ".venv/Scripts/python.exe"
SCRIPT = ROOT / "scripts/ankle_fe/ankle_lig_feb.py"
WORK = ROOT / "temp/thums/febio_lig"

BASE = {
    "DET_MIN": "0.0", "TRUNC_MM": "60", "AXIS_SRC": "full", "ALIGN_AXIS": "1",
    "FIX_TALUS": "all", "DRIVER": "dispx", "CART_NEW": "1", "TIE": "1",
    "CART_DIR": str(ROOT / "temp/thums/cart_m1"), "LOG_DATA": "1",
    "LIG_MODE": "none",
}


def parse(p):
    t = p.read_text(encoding="latin-1", errors="replace")
    cc = {int(m.group(1)): np.array([float(x) for x in m.group(2).split(",")])
          for m in re.finditer(r'<node id="(\d+)">([-\d.eE+,]+)</node>', t)}
    ns = {m.group(1): [int(x) for x in m.group(2).split(",") if x]
          for m in re.finditer(r'<NodeSet name="([^"]+)">([\d,]*)</NodeSet>', t)}
    return cc, ns


def last_disp(p):
    rows, tv, best = [], None, None
    for l in p.read_text(encoding="utf-8", errors="replace").split("\n"):
        if l.startswith("*Step"):
            if rows and tv is not None and (best is None or abs(tv - 1) < abs(best[0] - 1)):
                best = (tv, list(rows))
            rows, tv = [], None
        elif l.startswith("*Time"):
            tv = float(l.split("=")[1])
        elif l.strip() and not l.startswith("*"):
            try:
                rows.append([float(x) for x in l.split(",")])
            except ValueError:
                pass
    if rows and tv is not None and (best is None or abs(tv - 1) < abs(best[0] - 1)):
        best = (tv, rows)
    if best is None:
        return None
    return {int(r[0]): np.array(r[1:4]) for r in best[1]}


print(f"{'配置':<22} {'rc':>3} {'|反力z|(N)':>11} {'顶面dz(mm)':>11} "
      f"{'距骨dz':>9} {'相对接近':>10}  判定")
print("─" * 100)
for lbl, d in (("DISP=+1.0", "+1.0"), ("DISP=-1.0", "-1.0")):
    env = os.environ.copy()
    env.update(BASE)
    env["DISP_MM"] = d
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([str(PY), "-u", str(SCRIPT)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=3600)
    out = (r.stdout or "") + (r.stderr or "")
    rc = re.search(r"FEBio rc=(-?\d+)", out)
    cc, ns = parse(WORK / "ankle_lig_WK.feb")
    disp = last_disp(WORK / "disp_WK.txt") or {}
    tt = [i for i in ns["tibia_top_ns"] if i in disp and i in cc]
    tl = [i for i in ns["talus_low_ns"] if i in disp and i in cc]
    dz_t = float(np.mean([disp[i][2] for i in tt])) if tt else float("nan")
    dz_l = float(np.mean([disp[i][2] for i in tl])) if tl else float("nan")
    z_tib = float(np.mean([cc[i][2] for i in tt]))
    z_tal = float(np.mean([cc[i][2] for i in tl]))
    # 相对接近 = 顶面相对距骨沿 (距骨→顶面) 方向的位移；负 = 靠近 = 压实
    sgn = 1.0 if z_tib > z_tal else -1.0        # 顶面在距骨上方 ⇒ 接近 = dz 变负
    approach = sgn * (dz_t - dz_l)
    Fz = float("nan")
    rp = WORK / "reac_tibia_WK.txt"
    if rp.exists():
        rr = last_disp(rp)
        if rr:
            Fz = float(abs(sum(v[2] for v in rr.values())))
    verdict = ("**压实**（相对接近 ↑）" if approach > 1e-4 else
               "拉开（相对远离）" if approach < -1e-4 else "无相对运动")
    print(f"{lbl:<22} {rc.group(1) if rc else '?':>3} {Fz:>11.3f} {dz_t:>11.4f} "
          f"{dz_l:>9.4f} {approach:>10.4f}  {verdict}")
print(f"\n（参考：顶面 z̄={z_tib:.2f}，距骨 z̄={z_tal:.2f} ⇒ 距骨在"
      f"{'下' if z_tal < z_tib else '上'}方）")
