"""修正符号后的多步：落位=压缩 + 载荷=压缩 ⇒ 是否打通力控制。"""
import os
import re
import subprocess
from pathlib import Path

R = Path(r"D:\Project\climbing_fall_analysis")
W = R / "temp" / "thums" / "febio_lig"
FEB = r"D:\Program\FEBioStudio\bin\febio4.exe"
PY = str(R / ".venv" / "Scripts" / "python.exe")
CASES = [
    ("多步 压缩 不HOLD AUGLAG", dict(MULTISTEP="1", SEAT_MM="0.3",
                                  SEAT_HOLD="0", LAUGON="AUGLAG")),
    ("多步 压缩 不HOLD PENALTY", dict(MULTISTEP="1", SEAT_MM="0.3",
                                    SEAT_HOLD="0", LAUGON="PENALTY")),
    ("单步 压缩 AUGLAG", dict(LAUGON="AUGLAG")),
    ("单步 压缩 PENALTY", dict(LAUGON="PENALTY")),
]
for tag, ov in CASES:
    e = dict(os.environ)
    e.update(DRIVER="force", FORCE_N="600", TRUNC_MM="60", AXIS_SRC="full",
             ALIGN_AXIS="1", CART_DIR=str(R / "temp" / "thums" / "cart_m1"),
             CART_NEW="1", TIE="1", DET_MIN="0.0", TALUS_RIGID="1",
             TAL_ROTLOCK="0", FIX_TALUS="transverse", TIME_STEPS="20",
             LOAD_SIGN="+1", TIME_STEPS_SEAT="5")
    e.update(ov)
    p = subprocess.run([PY, "-u",
                        str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
                       capture_output=True, text=True, env=e, cwd=str(R),
                       timeout=900)
    if p.returncode:
        print(f"{tag:<26} build rc={p.returncode}", flush=True)
        continue
    r = subprocess.run([FEB, "-i", "ankle_lig_WK.feb", "-o", "s.xplt"],
                       capture_output=True, text=True, cwd=str(W), timeout=1800)
    out = (r.stdout or "") + (r.stderr or "")
    fin = [float(x) for x in re.findall(r"converged at time : ([\d.]+)", out)]
    mx = max(fin) if fin else None
    print(f"{tag:<26} t={mx if mx else '-':<9} 载荷≈"
          f"{round(600*(mx-1),1) if mx and mx > 1 else (round(600*mx,1) if mx else 0):>6} N"
          f"  {'NORMAL ✅' if 'N O R M A L' in out else 'FAIL ❌'}", flush=True)
