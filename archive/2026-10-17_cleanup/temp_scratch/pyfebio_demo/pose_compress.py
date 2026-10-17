"""压缩侧姿势矩阵（DISP_MM=-1.0 ⇒ 约 -609 N 轴向），输出三分量合力。"""
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
print(f"{'TILT':>5} {'TILT_ML':>8} {'NORMAL':>7} {'Fx':>9} {'Fy':>9} {'Fz':>9}"
      f" {'|Fshear|':>9}")
for a, b in CASES:
    e = dict(os.environ)
    e.update(DRIVER="dispx", DISP_MM="-1.0", TRUNC_MM="60", AXIS_SRC="full",
             ALIGN_AXIS="1", CART_DIR=str(R / "temp" / "thums" / "cart_m1"),
             CART_NEW="1", TIE="1", DET_MIN="0.0", TALUS_RIGID="1",
             TAL_ROTLOCK="0", FIX_TALUS="transverse", TIME_STEPS="20",
             LOG_DATA="1", LAUGON="PENALTY", POSE_TILT=str(a),
             POSE_TILT_ML=str(b))
    subprocess.run([PY, "-u",
                    str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
                   capture_output=True, text=True, env=e, cwd=str(R), timeout=900)
    r = subprocess.run([FEB, "-i", "ankle_lig_WK.feb", "-o", "s.xplt"],
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
    print(f"{a:>5} {b:>8} {str(ok):>7} {G[0]:>9.2f} {G[1]:>9.2f} {G[2]:>9.2f}"
          f" {np.hypot(G[0],G[1]):>9.2f}", flush=True)
