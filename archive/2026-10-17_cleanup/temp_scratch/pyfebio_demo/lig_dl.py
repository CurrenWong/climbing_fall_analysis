"""逐韧带 ΔL：真实内外翻姿势下的韧带伸长（判据 C3/C4 的核心量）。"""
import io
import os
import re
import subprocess
from pathlib import Path

import numpy as np

R = Path(r"D:\Project\climbing_fall_analysis")
W = R / "temp" / "thums" / "febio_lig"
FEB = r"D:\Program\FEBioStudio\bin\febio4.exe"
PY = str(R / ".venv" / "Scripts" / "python.exe")
COMMON = dict(DRIVER="dispx", DISP_MM="-1.0", LOAD_SIGN="+1",
              FIX_TALUS="transverse", TRUNC_MM="60", AXIS_SRC="full",
              ALIGN_AXIS="1", CART_NEW="1", TIE="1", DET_MIN="0.0",
              TALUS_RIGID="0", TAL_ROTLOCK="0", TIME_STEPS="20",
              LOG_DATA="1", LAUGON="PENALTY",
              CART_DIR=str(R / "temp" / "thums" / "cart_m1"))


def coords_of(fp):
    t = fp.read_text(encoding="latin-1", errors="replace")
    X = {}
    for m in re.finditer(r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)'
                         r'\s*,\s*([-\d.eE+]+)', t):
        X[int(m.group(1))] = [float(m.group(2)), float(m.group(3)),
                              float(m.group(4))]
    return X


def disp_last(fp):
    rows, tmax = [], []
    for ln in io.open(fp, encoding="utf-8", errors="replace"):
        ln = ln.strip()
        if ln.startswith("*Step"):
            tmax.append([])
            continue
        if not ln or ln.startswith("*"):
            continue
        try:
            v = [float(x) for x in ln.replace(",", " ").split()]
        except ValueError:
            continue
        if tmax:
            tmax[-1].append(v)
    last = tmax[-1] if tmax else []
    return {int(v[0]): np.array(v[1:4]) for v in last}


res = {}
for ie in ("-3", "-1", "0", "+1", "+3", "+5"):
    feb = W / "ankle_lig_WK.feb"
    if feb.exists():
        feb.unlink()
    e = dict(os.environ)
    e.update(COMMON)
    e["POSE_IE"] = ie
    subprocess.run([PY, "-u", str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
                   capture_output=True, text=True, env=e, cwd=str(R), timeout=900)
    if not feb.exists():
        print(f"IE={ie}: 构建失败")
        continue
    subprocess.run([FEB, "-i", feb.name, "-o", "s.xplt"], capture_output=True,
                   text=True, cwd=str(W), timeout=1800)
    X = coords_of(feb)
    U = disp_last(W / "disp_WK.txt")
    z = np.load(W / "lig_read_WK.npz", allow_pickle=True)
    names = [str(x) for x in z["lig_names"]]
    PR = z["lig_pairs"]
    npr = z["lig_npair"]
    out = {}
    for i, nm in enumerate(names):
        k = int(npr[i])
        sub = PR[i, :k, :]                      # (10, 13, 2) 布局；-1 为填充
        dl = []
        for a, b in sub:
            a, b = int(a), int(b)
            if a < 0 or b < 0 or a not in X or b not in X:
                continue
            p0 = np.array(X[a]) - np.array(X[b])
            ua = U.get(a, np.zeros(3))
            ub = U.get(b, np.zeros(3))
            p1 = p0 + (ua - ub)
            L0 = np.linalg.norm(p0)
            if L0 < 1e-9:
                continue
            dl.append((np.linalg.norm(p1) - L0) / L0 * 100.0)
        out[nm] = (float(np.mean(dl)) if dl else float("nan"),
                   float(np.max(dl)) if dl else float("nan"), len(dl))
    res[ie] = out
    print(f"\n=== POSE_IE={ie} ===（应变 % ：均值 / 最大 / 对数）")
    for nm, (mu, mx, n) in sorted(out.items(), key=lambda x: -abs(x[1][1] or 0)):
        print(f"  {nm:<28} {mu:>8.3f} {mx:>8.3f}  n={n}")

print("\n===== 汇总：各姿势的最大韧带应变 (%) =====")
allnm = sorted({nm for o in res.values() for nm in o})
hdr = "  ".join(f"{ie:>7}" for ie in res)
print(f"{'ligament':<28} {hdr}")
for nm in allnm:
    vals = "  ".join(f"{res[ie].get(nm, (float('nan'),))[1]:>7.3f}" for ie in res)
    print(f"{nm:<28} {vals}")
