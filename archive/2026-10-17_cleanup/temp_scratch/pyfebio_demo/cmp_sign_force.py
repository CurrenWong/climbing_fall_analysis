"""决定性：HOLD+AUGLAG 下对比两个符号的末步反力（应 ≈ 座力+600 = 665 N）。"""
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


for sgn in ("-1", "+1"):
    e = dict(os.environ)
    e.update(DRIVER="force", FORCE_N="600", TRUNC_MM="60", AXIS_SRC="full",
             ALIGN_AXIS="1", CART_DIR=str(R / "temp" / "thums" / "cart_m1"),
             CART_NEW="1", TIE="1", DET_MIN="0.0", TALUS_RIGID="1",
             TAL_ROTLOCK="0", FIX_TALUS="transverse", TIME_STEPS="20",
             LAUGON="AUGLAG", MULTISTEP="1", SEAT_MM="0.1", SEAT_HOLD="1",
             TIME_STEPS_SEAT="5", LOAD_SIGN=sgn, LOG_DATA="1")
    subprocess.run([PY, "-u",
                    str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
                   capture_output=True, text=True, env=e, cwd=str(R), timeout=900)
    r = subprocess.run([FEB, "-i", "ankle_lig_WK.feb", "-o", "ms.xplt"],
                       capture_output=True, text=True, cwd=str(W), timeout=1800)
    out = (r.stdout or "") + (r.stderr or "")
    G = np.zeros(3)
    for nm in ("tibia", "fibula"):
        a = rows(next(W.glob(f"reac_{nm}_*.txt")))
        n = len(np.unique(a[:, 0].astype(int)))
        F = a[:, 1:4].reshape(len(a) // n, n, 3).sum(axis=1)
        G = G + F[-1]
    print(f"LOAD_SIGN={sgn:>2}  NORMAL={'N O R M A L' in out}  "
          f"末步顶面合力 {np.round(G,2)}  |F|={np.linalg.norm(G):.1f} N", flush=True)
