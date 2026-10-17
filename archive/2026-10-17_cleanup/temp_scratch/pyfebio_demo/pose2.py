"""压缩侧姿势矩阵（修好 _outer 后；带 rc 检查 + 删旧 .feb 防读到脏文件）。"""
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
FEBF = W / "ankle_lig_WK.feb"


def rows(fp):
    o = []
    for ln in io.open(fp, encoding="utf-8", errors="replace"):
        ln = ln.strip()
        if ln and not ln.startswith("*"):
            try:
                o.append([float(x) for x in ln.replace(",", " ").split()])
            except ValueError:
                pass
    return np.array(o)


CASES = [(0, 0), (10, 0), (-10, 0), (20, 0), (0, 10), (0, -10), (10, 10), (10, -10)]
print(f"{'TILT':>5} {'T_ML':>5} {'NORMAL':>7} {'Fx':>9} {'Fy':>9} {'Fz':>9}"
      f" {'Fshear':>8}", flush=True)
for a, b in CASES:
    if FEBF.exists():
        FEBF.unlink()
    e = dict(os.environ)
    e.update(DRIVER="dispx", DISP_MM="-1.0", TRUNC_MM="60", AXIS_SRC="full",
             ALIGN_AXIS="1", CART_DIR=str(R / "temp" / "thums" / "cart_m1"),
             CART_NEW="1", TIE="1", DET_MIN="0.0", TALUS_RIGID="1",
             TAL_ROTLOCK="0", FIX_TALUS="transverse", TIME_STEPS="20",
             LOG_DATA="1", LAUGON="PENALTY", POSE_TILT=str(a),
             POSE_TILT_ML=str(b))
    p = subprocess.run([PY, "-u",
                        str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
                       capture_output=True, text=True, env=e, cwd=str(R), timeout=900)
    if p.returncode or not FEBF.exists():
        print(f"{a:>5} {b:>5}   build rc={p.returncode}", flush=True)
        continue
    r = subprocess.run([FEB, "-i", FEBF.name, "-o", "s.xplt"],
                       capture_output=True, text=True, cwd=str(W), timeout=1800)
    out = (r.stdout or "") + (r.stderr or "")
    G = np.zeros(3)
    for nm in ("tibia", "fibula"):
        try:
            arr = rows(next(W.glob(f"reac_{nm}_*.txt")))
        except StopIteration:
            continue
        n = len(np.unique(arr[:, 0].astype(int)))
        G = G + arr[:, 1:4].reshape(len(arr) // n, n, 3).sum(axis=1)[-1]
    ok = "N O R M A L" in out
    print(f"{a:>5} {b:>5} {str(ok):>7} {G[0]:>9.2f} {G[1]:>9.2f} {G[2]:>9.2f}"
          f" {np.hypot(G[0],G[1]):>8.2f}", flush=True)
