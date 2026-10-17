"""用可信方法验收 cart_m1_fix：建 .feb → 在 .feb 上量 worst-det(J)。"""
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

R = Path(r"D:\Project\climbing_fall_analysis")
GP2 = 1 / np.sqrt(3)


def run(cart_dir):
    env = dict(os.environ)
    env.update({"DRIVER": "dispx", "DISP_MM": "+1.6", "DET_MIN": "0.0",
                "TRUNC_MM": "60", "AXIS_SRC": "full", "ALIGN_AXIS": "1",
                "FIX_TALUS": "transverse", "CART_NEW": "1", "TIE": "1",
                "CART_DIR": cart_dir, "POSE_IE": "0",
                "PYTHONIOENCODING": "utf-8"})
    r = subprocess.run([str(R / ".venv/Scripts/python.exe"), "-u",
                        str(R / "scripts/ankle_fe/ankle_lig_feb.py")],
                       cwd=str(R), env=env, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=3600)
    o = (r.stdout or "") + (r.stderr or "")
    return o


def measure(febpath):
    txt = Path(febpath).read_text(encoding="latin-1", errors="replace")
    XYZ = {}
    for m in re.finditer(r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,'
                         r'\s*([-\d.eE+]+)', txt):
        XYZ[int(m.group(1))] = np.array([float(m.group(2)), float(m.group(3)),
                                         float(m.group(4))])
    DOMS = {}
    for m in re.finditer(r'<Elements(\s[^>]*?)>(.*?)</Elements>', txt, re.S):
        a, body = m.group(1), m.group(2)
        k = re.search(r'type="(\w+)"', a)
        n = re.search(r'name="([^"]*)"', a)
        if not (k and n):
            continue
        DOMS[n.group(1)] = (k.group(1), [[int(x) for x in c.split(",")] for c in
                                         re.findall(r'<elem[^>]*>\s*([\d,\s]+?)\s*'
                                                    r'</elem>', body)])

    def dNw(r, s, t):
        return np.array([-(1-t), -(1-t), -(1-r-s), (1-t), 0, -r, 0, (1-t), -s,
                         -t, -t, (1-r-s), t, 0, r, 0, t, s]).reshape(6, 3)

    RS = (0.0, 1/6, 1/3, 1/2, 2/3, 1.0)
    TS = (0.0, 0.25, 0.5, 0.75, 1.0)
    out = {}
    for nm, (kind, conn) in DOMS.items():
        if kind != "penta6" or "tibiotalar" not in nm:
            continue
        ws, th = [], []
        for c in conn:
            P = np.array([XYZ[n] for n in c], float)
            w = min(float(np.linalg.det(dNw(r, s, t).T @ P))
                    for r in RS for s in RS if r + s <= 1 + 1e-9 for t in TS)
            ws.append(w)
            th.append(float(np.linalg.norm(P[[3, 4, 5]].mean(0)
                                           - P[[0, 1, 2]].mean(0))))
        ws, th = np.array(ws), np.array(th)
        out[nm] = (len(ws), (ws <= 0).sum(), (ws <= 1e-2).sum(), ws.min(),
                   th.min(), np.median(th), th.max(), (th > 1.6).sum())
    return out


RES = {}
for label, cd in (("基线 cart_m1", "temp/thums/cart_m1"),
                  ("新 cart_m1_fix", "temp/thums/cart_m1_fix")):
    print("=" * 100)
    print(f"=== {label}")
    print("=" * 100)
    o = run(cd)
    for l in o.splitlines():
        if re.search(r"FEBio rc|外表面顶|area自检·完整|保留.*penta6|截断", l):
            print("  " + l.strip()[:140])
    m = re.search(r"两顶面合力.*?\|F\|=([\d.]+)", o)
    F = float(m.group(1)) if m else float("nan")
    print(f"  → 合力 = {F:.4f} N")
    # ★ 立刻量这一次的 .feb（否则被下一次覆盖）
    feb = R / "temp/thums/febio_lig/ankle_lig_WK.feb"
    d = measure(feb)
    for nm, (n, n0, n2, wmin, tmin, tmed, tmax, nov) in sorted(d.items()):
        print(f"      {nm:<30} 单元 {n:>4}  ≤0 {n0:>4}  ≤1e-2 {n2:>4}  "
              f"worst-min {wmin:>12.4e}")
        print(f"          厚度 min {tmin:.4f} 中位 {tmed:.4f} max {tmax:.4f} mm "
              f"(>1.6 的 {nov} 个)")
    RES[label] = (F, d)

print("\n" + "=" * 100)
print("对照汇总")
print("=" * 100)
for label, (F, d) in RES.items():
    tot0 = sum(v[1] for v in d.values())
    tot2 = sum(v[2] for v in d.values())
    nov = sum(v[7] for v in d.values())
    print(f"  {label:<18} F={F:>9.4f} N  ≤0 合计 {tot0:>4}  ≤1e-2 合计 {tot2:>4}"
          f"  厚度>1.6 的 {nov:>4}")
