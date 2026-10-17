"""独立裁决：真·压实方向 ⇒ 关节面间距必须**变小**。

不看任何符号约定，只看几何：取两个软骨域的节点，算 t=0 / t=1 的最小间距。
"""
import io
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


def doms(p):
    """从 .feb 取 tibiotalar_tibia / tibiotalar_talus 两个域的节点。"""
    t = p.read_text(encoding="latin-1", errors="replace")
    out = {}
    for m in re.finditer(
            r'<SolidDomain name="(tibiotalar_\w+)__\w+"[^>]*>(.*?)</SolidDomain>',
            t, re.S):
        ids = sorted({int(x) for x in re.findall(r"\d+", m.group(2))})
        out[m.group(1)] = ids
    return out


def coords(p):
    t = p.read_text(encoding="latin-1", errors="replace")
    return {int(m.group(1)): np.array([float(x) for x in m.group(2).split(",")])
            for m in re.finditer(r'<node id="(\d+)">([-\d.eE+,]+)</node>', t)}


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


print(f"{'配置':<28} {'rc':>3} {'|F|(N)':>10} {'间距 t=0':>10} {'间距 t=1':>10} "
      f"{'Δ':>9}  判定")
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
    feb = WORK / "ankle_lig_WK.feb"
    dm = doms(feb)
    cc = coords(feb)
    disp = last_disp(WORK / "disp_WK.txt")
    A = np.array([cc[i] for i in dm["tibiotalar_tibia"]])
    B = np.array([cc[i] for i in dm["tibiotalar_talus"]])
    d0 = float(np.sqrt(((A[:, None, :] - B[None, :, :]) ** 2).sum(-1)).min())
    A1 = np.array([cc[i] + disp[i] for i in dm["tibiotalar_tibia"]] if disp else A)
    B1 = np.array([cc[i] + disp[i] for i in dm["tibiotalar_talus"]] if disp else B)
    d1 = float(np.sqrt(((A1[:, None, :] - B1[None, :, :]) ** 2).sum(-1)).min())
    f = re.search(r"FEBio rc=(-?\d+)", out)
    # 用反力文件算合力
    F = float("nan")
    rp = WORK / "reac_tibia_WK.txt"
    if rp.exists():
        rr = last_disp(rp)
        if rr:
            F = float(abs(sum(v[2] for v in rr.values())))
    verdict = "**压实**（间距↓）" if d1 < d0 - 1e-4 else (
        "拉开（间距↑）" if d1 > d0 + 1e-4 else "无变化")
    print(f"{lbl:<28} {f.group(1) if f else '?':>3} {F:>10.3f} "
          f"{d0:>10.4f} {d1:>10.4f} {d1-d0:>9.4f}  {verdict}")
print("\n★ 这是**独立于任何符号约定**的几何裁决。")
