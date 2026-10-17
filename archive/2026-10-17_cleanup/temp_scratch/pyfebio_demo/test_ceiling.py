"""12 N 天花板：是哪个构件先顶不住？"""
import os
import re
import subprocess
from pathlib import Path

R = Path(r"D:\Project\climbing_fall_analysis")
W = R / "temp" / "thums" / "febio_lig"
FEB = r"D:\Program\FEBioStudio\bin\febio4.exe"
CASES = [
    ("基线 600N 胫+腓", dict(FORCE_N="600")),
    ("腓骨不分担", dict(FORCE_N="600", FIB_SHARE="0")),
    ("100N 量级", dict(FORCE_N="100")),
    ("5000N 量级", dict(FORCE_N="5000")),
    ("无韧带 LIG_MODE=none", dict(FORCE_N="600", LIG_MODE="none")),
]
for tag, ov in CASES:
    e = dict(os.environ)
    e.update(DRIVER="force", FIX_TALUS="transverse", TRUNC_MM="60",
             AXIS_SRC="full", ALIGN_AXIS="1",
             CART_DIR=str(R / "temp" / "thums" / "cart_m1"), CART_NEW="1",
             TIE="1", DET_MIN="0.0", TALUS_RIGID="1", TAL_ROTLOCK="0",
             MULTISTEP="1", SEAT_MM="0.10", TIME_STEPS="40", TIME_STEPS_SEAT="20")
    e.update(ov)
    subprocess.run([str(R / ".venv" / "Scripts" / "python.exe"), "-u",
                    str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
                   capture_output=True, text=True, env=e, cwd=str(R), timeout=3000)
    r = subprocess.run([FEB, "-i", "ankle_lig_WK.feb", "-o", "ms.xplt"],
                       capture_output=True, text=True, cwd=str(W), timeout=3000)
    out = (r.stdout or "") + (r.stderr or "")
    fin = [float(x) for x in re.findall(r"converged at time : ([\d.]+)", out)]
    fail = re.findall(r"failed to converge at time : ([\d.]+)", out)
    fN = float(ov.get("FORCE_N", "600"))
    s2 = [t for t in fin if t > 1.0]
    L = round(fN * (max(s2) - 1), 1) if s2 else 0
    print(f"{tag:<20} step2 到 t={max(s2) if s2 else '-':<10} 载荷 {L:>7.1f} N / "
          f"{fN:>6.0f} N  {('NORMAL ✅' if 'N O R M A L' in out else 'FAIL ❌')}")
