"""关键符号判据：dispx 592.61 N 是压（-Rz）还是拉（+Rz）？并测反向位移。"""
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


for tag, ov in (("dispx +1.6 (既定)", dict(DRIVER="dispx", DISP_MM="1.6")),
                ("dispx -1.6 (反向)", dict(DRIVER="dispx", DISP_MM="-1.6"))):
    e = dict(os.environ)
    e.update(FIX_TALUS="transverse", TRUNC_MM="60", AXIS_SRC="full",
             ALIGN_AXIS="1", CART_DIR=str(R / "temp" / "thums" / "cart_m1"),
             CART_NEW="1", TIE="1", DET_MIN="0.0", TALUS_RIGID="1",
             TAL_ROTLOCK="0", TIME_STEPS="20", LOG_DATA="1")
    e.update(ov)
    p = subprocess.run([PY, "-u",
                        str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
                       capture_output=True, text=True, env=e, cwd=str(R),
                       timeout=900)
    if p.returncode:
        print(f"{tag:<20} build rc={p.returncode}", flush=True)
        continue
    r = subprocess.run([FEB, "-i", "ankle_lig_WK.feb", "-o", "s.xplt"],
                       capture_output=True, text=True, cwd=str(W), timeout=1800)
    out = (r.stdout or "") + (r.stderr or "")
    G = np.zeros(3)
    for nm in ("tibia", "fibula"):
        a = rows(next(W.glob(f"reac_{nm}_*.txt")))
        n = len(np.unique(a[:, 0].astype(int)))
        F = a[:, 1:4].reshape(len(a) // n, n, 3).sum(axis=1)
        G = G + F[-1]
    print(f"{tag:<20} NORMAL={'N O R M A L' in out}  末步顶面合力 {np.round(G,2)}"
          f"  Fz={G[2]:+.2f} N  ⇒ {'压缩(-z)' if G[2] < 0 else '拉伸(+z)'}"
          if G[2] != 0 else f"{tag}: Fz=0", flush=True)
