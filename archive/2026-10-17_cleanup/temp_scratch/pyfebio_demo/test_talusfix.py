"""距骨约束方式 × 多步：哪个组合能真正把力加上去。"""
import os
import re
import subprocess
from pathlib import Path

R = Path(r"D:\Project\climbing_fall_analysis")
W = R / "temp" / "thums" / "febio_lig"
FEB = r"D:\Program\FEBioStudio\bin\febio4.exe"
CASES = [
    ("transverse + HOLD", dict(FIX_TALUS="transverse", SEAT_HOLD="1")),
    ("transverse 不 HOLD", dict(FIX_TALUS="transverse", SEAT_HOLD="0")),
    ("all + HOLD", dict(FIX_TALUS="all", SEAT_HOLD="1")),
    ("none + HOLD", dict(FIX_TALUS="none", SEAT_HOLD="1")),
    ("transverse+ROTLOCK + HOLD", dict(FIX_TALUS="transverse", SEAT_HOLD="1",
                                      TAL_ROTLOCK="1")),
]
print(f"{'case':<28} {'step2 到':>9} {'载荷 N':>8}  结果")
for tag, ov in CASES:
    e = dict(os.environ)
    e.update(DRIVER="force", FORCE_N="600", TRUNC_MM="60", AXIS_SRC="full",
             ALIGN_AXIS="1", CART_DIR=str(R / "temp" / "thums" / "cart_m1"),
             CART_NEW="1", TIE="1", DET_MIN="0.0", TALUS_RIGID="1",
             TAL_ROTLOCK="0", MULTISTEP="1", SEAT_MM="0.10",
             TIME_STEPS="40", TIME_STEPS_SEAT="20")
    e.update(ov)
    subprocess.run([str(R / ".venv" / "Scripts" / "python.exe"), "-u",
                    str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
                   capture_output=True, text=True, env=e, cwd=str(R), timeout=3000)
    r = subprocess.run([FEB, "-i", "ankle_lig_WK.feb", "-o", "ms.xplt"],
                       capture_output=True, text=True, cwd=str(W), timeout=3000)
    out = (r.stdout or "") + (r.stderr or "")
    fin = [float(x) for x in re.findall(r"converged at time : ([\d.]+)", out)]
    s2 = [t for t in fin if t > 1.0]
    mx = max(s2) if s2 else None
    L = round(600 * (mx - 1), 1) if mx else 0.0
    print(f"{tag:<28} {mx if mx else '-':>9} {L:>8.1f}  "
          f"{'NORMAL ✅' if 'N O R M A L' in out else 'FAIL ❌'}")
