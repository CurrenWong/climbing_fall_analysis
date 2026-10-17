"""扫小落位量：找"刚好接触"的甜点（step1 末反力小 ⇒ 撤约束时回弹小）。"""
import io
import os
import subprocess
from pathlib import Path

import numpy as np

R = Path(r"D:\Project\climbing_fall_analysis")
W = R / "temp" / "thums" / "febio_lig"
FEB = r"D:\Program\FEBioStudio\bin\febio4.exe"


def run(seat_pct, ts):
    e = dict(os.environ)
    e.update(DRIVER="force", FORCE_N="600", FIX_TALUS="transverse", TRUNC_MM="60",
             AXIS_SRC="full", ALIGN_AXIS="1",
             CART_DIR=str(R / "temp" / "thums" / "cart_m1"), CART_NEW="1",
             TIE="1", DET_MIN="0.0", TALUS_RIGID="1", TAL_ROTLOCK="0",
             MULTISTEP="1", SEAT_MM=str(seat_pct), TIME_STEPS=str(ts),
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
    seat_end = F[min(19, len(F) - 1), 2]
    norm = "N O R M A L" in out
    last = F[-1, 2]
    return norm, len(F), seat_end, last


print(f"{'SEAT_MM':>8} {'步数':>5} {'step1末Rz':>11} {'末步Rz':>10}  结果")
for seat in ("0.02", "0.05", "0.10", "0.20", "0.30", "0.50"):
    norm, ns, s1, last = run(seat, 20)
    print(f"{seat:>8} {ns:>5} {s1:>11.2f} {last:>10.2f}  "
          f"{'NORMAL ✅' if norm else 'FAIL ❌'}")
