"""THE 关键测试：不保持约束（载荷接管）+ AUGLAG ⇒ 反力应 ≈600 N。"""
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


def run(ov, tag):
    e = dict(os.environ)
    e.update(DRIVER="force", FORCE_N="600", TRUNC_MM="60", AXIS_SRC="full",
             ALIGN_AXIS="1", CART_DIR=str(R / "temp" / "thums" / "cart_m1"),
             CART_NEW="1", TIE="1", DET_MIN="0.0", TALUS_RIGID="1",
             TAL_ROTLOCK="0", MULTISTEP="1", SEAT_MM="0.10",
             FIX_TALUS="transverse", TIME_STEPS="20", TIME_STEPS_SEAT="5",
             LOG_DATA="1")
    e.update(ov)
    subprocess.run([PY, "-u",
                    str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
                   capture_output=True, text=True, env=e, cwd=str(R), timeout=900)
    r = subprocess.run([FEB, "-i", "ankle_lig_WK.feb", "-o", "ms.xplt"],
                       capture_output=True, text=True, cwd=str(W), timeout=1800)
    out = (r.stdout or "") + (r.stderr or "")

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

    G = np.zeros(3)
    for nm in ("tibia", "fibula"):
        a = rows(next(W.glob(f"reac_{nm}_*.txt")))
        n = len(np.unique(a[:, 0].astype(int)))
        F = a[:, 1:4].reshape(len(a) // n, n, 3).sum(axis=1)
        G = G + F[-1]
    fin = [float(x) for x in re.findall(r"converged at time : ([\d.]+)", out)]
    print(f"{tag:<28} NORMAL={'N O R M A L' in out}  末收敛 t={max(fin) if fin else '-'}"
          f"  末步顶面合力={np.round(G,2)} |F|={np.linalg.norm(G):.2f} N", flush=True)


run(dict(SEAT_HOLD="0", LAUGON="AUGLAG"), "不HOLD + AUGLAG")
run(dict(SEAT_HOLD="1", LAUGON="AUGLAG"), "HOLD + AUGLAG (对照)")
