"""决定性判据：翻转载荷符号（LOAD_SIGN=+1）后力控制是否通。"""
import os
import re
import subprocess
from pathlib import Path

R = Path(r"D:\Project\climbing_fall_analysis")
W = R / "temp" / "thums" / "febio_lig"
FEB = r"D:\Program\FEBioStudio\bin\febio4.exe"
PY = str(R / ".venv" / "Scripts" / "python.exe")


def go(tag, ov):
    e = dict(os.environ)
    e.update(DRIVER="force", FORCE_N="600", TRUNC_MM="60", AXIS_SRC="full",
             ALIGN_AXIS="1", CART_DIR=str(R / "temp" / "thums" / "cart_m1"),
             CART_NEW="1", TIE="1", DET_MIN="0.0", TALUS_RIGID="1",
             TAL_ROTLOCK="0", FIX_TALUS="transverse", TIME_STEPS="20",
             LAUGON="AUGLAG")
    e.update(ov)
    p = subprocess.run([PY, "-u",
                        str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
                       capture_output=True, text=True, env=e, cwd=str(R),
                       timeout=900)
    if p.returncode:
        print(f"{tag:<30} build rc={p.returncode}", flush=True)
        return
    r = subprocess.run([FEB, "-i", "ankle_lig_WK.feb", "-o", "s.xplt"],
                       capture_output=True, text=True, cwd=str(W), timeout=1800)
    out = (r.stdout or "") + (r.stderr or "")
    fin = [float(x) for x in re.findall(r"converged at time : ([\d.]+)", out)]
    mx = max(fin) if fin else None
    print(f"{tag:<30} t={mx if mx else '-':<9} "
          f"{'NORMAL ✅' if 'N O R M A L' in out else 'FAIL ❌'}", flush=True)


for sgn in ("-1", "+1"):
    go(f"单步 LOAD_SIGN={sgn}  AUGLAG", dict(LOAD_SIGN=sgn))
for sgn in ("-1", "+1"):
    go(f"多步 LOAD_SIGN={sgn}  AUGLAG", dict(LOAD_SIGN=sgn, MULTISTEP="1",
                                           SEAT_MM="0.1", SEAT_HOLD="1",
                                           TIME_STEPS_SEAT="5"))
