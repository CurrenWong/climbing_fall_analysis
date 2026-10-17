"""对照：step2 几乎不加载时，Rz 是否仍保持 step1 的 245 N（= 位移仍被保持）。"""
import io
import os
import subprocess
from pathlib import Path

import numpy as np

R = Path(r"D:\Project\climbing_fall_analysis")
W = R / "temp" / "thums" / "febio_lig"
FEB = r"D:\Program\FEBioStudio\bin\febio4.exe"


def run(force, tag):
    e = dict(os.environ)
    e.update(DRIVER="force", FORCE_N=str(force), FIX_TALUS="transverse",
             TRUNC_MM="60", AXIS_SRC="full", ALIGN_AXIS="1",
             CART_DIR=str(R / "temp" / "thums" / "cart_m1"), CART_NEW="1",
             TIE="1", DET_MIN="0.0", TALUS_RIGID="1", TAL_ROTLOCK="0",
             MULTISTEP="1", SEAT_MM="0.5", TIME_STEPS="20",
             TIME_STEPS_SEAT="20", LOG_DATA="1")
    subprocess.run([str(R / ".venv" / "Scripts" / "python.exe"), "-u",
                    str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
                   capture_output=True, text=True, env=e, cwd=str(R), timeout=3000)
    r = subprocess.run([FEB, "-i", "ankle_lig_WK.feb", "-o", "ms.xplt"],
                       capture_output=True, text=True, cwd=str(W), timeout=3000)
    out = (r.stdout or "") + (r.stderr or "")
    rows = []
    for ln in io.open(W / "reac_tibia_WK.txt", encoding="utf-8", errors="replace"):
        ln = ln.strip()
        if ln and not ln.startswith("*"):
            try:
                rows.append([float(x) for x in ln.replace(",", " ").split()])
            except ValueError:
                pass
    a = np.array(rows)
    n = len(np.unique(a[:, 0].astype(int)))
    F = a[:, 1:4].reshape(len(a) // n, n, 3).sum(axis=1)
    norm = "N O R M A L" in out
    fin = [l for l in out.splitlines() if "converged at time" in l]
    print(f"{tag}: 时间点 {len(F)}  {'NORMAL ✅' if norm else 'FAIL ❌'}")
    print(f"   Rz: {np.round(F[:, 2], 2)}")


run(0.001, "FORCE_N≈0")
