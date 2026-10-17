"""干净的双向符号判据（小位移，都应收敛）。"""
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


for d in ("-0.25", "+0.25"):
    e = dict(os.environ)
    e.update(DRIVER="dispx", DISP_MM=d, TRUNC_MM="60", AXIS_SRC="full",
             ALIGN_AXIS="1", CART_DIR=str(R / "temp" / "thums" / "cart_m1"),
             CART_NEW="1", TIE="1", DET_MIN="0.0", TALUS_RIGID="1",
             TAL_ROTLOCK="0", FIX_TALUS="transverse", TIME_STEPS="20",
             LOG_DATA="1", LAUGON="PENALTY")
    subprocess.run([PY, "-u",
                    str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
                   capture_output=True, text=True, env=e, cwd=str(R), timeout=900)
    r = subprocess.run([FEB, "-i", "ankle_lig_WK.feb", "-o", "s.xplt"],
                       capture_output=True, text=True, cwd=str(W), timeout=1800)
    out = (r.stdout or "") + (r.stderr or "")
    G = np.zeros(3)
    for nm in ("tibia", "fibula"):
        a = rows(next(W.glob(f"reac_{nm}_*.txt")))
        n = len(np.unique(a[:, 0].astype(int)))
        G = G + a[:, 1:4].reshape(len(a) // n, n, 3).sum(axis=1)[-1]
    print(f"DISP_MM={d:>6}  NORMAL={'N O R M A L' in out}  "
          f"Fz={G[2]:+9.2f} N  ⇒ {'压缩 ✅' if G[2] < 0 else '拉伸 ❌'}", flush=True)
