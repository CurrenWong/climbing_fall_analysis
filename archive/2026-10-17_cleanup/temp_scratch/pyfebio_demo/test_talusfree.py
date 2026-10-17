"""快测：距骨轴向也自由（只靠弹簧）⇒ 两端都弹性支撑，避免过约束。"""
import os
import re
import subprocess
from pathlib import Path

R = Path(r"D:\Project\climbing_fall_analysis")
W = R / "temp" / "thums" / "febio_lig"
FEB = r"D:\Program\FEBioStudio\bin\febio4.exe"
PY = str(R / ".venv" / "Scripts" / "python.exe")
CASES = [
    ("none k=500 单步", dict(FIX_TALUS="none", LAUGON="AUGLAG")),
    ("none k=5000 单步", dict(FIX_TALUS="none", LAUGON="AUGLAG",
                            TALUS_SPRING_K="5000")),
    ("none k=50000 单步", dict(FIX_TALUS="none", LAUGON="AUGLAG",
                              TALUS_SPRING_K="50000")),
    ("none k=50000 PENALTY", dict(FIX_TALUS="none")),
]
for tag, ov in CASES:
    e = dict(os.environ)
    e.update(DRIVER="force", FORCE_N="600", TRUNC_MM="60", AXIS_SRC="full",
             ALIGN_AXIS="1", CART_DIR=str(R / "temp" / "thums" / "cart_m1"),
             CART_NEW="1", TIE="1", DET_MIN="0.0", TALUS_RIGID="1",
             TAL_ROTLOCK="0", TIME_STEPS="20")
    e.update(ov)
    p = subprocess.run([PY, "-u",
                        str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
                       capture_output=True, text=True, env=e, cwd=str(R),
                       timeout=900)
    if p.returncode:
        print(f"{tag:<24} build rc={p.returncode}", flush=True)
        continue
    r = subprocess.run([FEB, "-i", "ankle_lig_WK.feb", "-o", "s.xplt"],
                       capture_output=True, text=True, cwd=str(W), timeout=1800)
    out = (r.stdout or "") + (r.stderr or "")
    fin = [float(x) for x in re.findall(r"converged at time : ([\d.]+)", out)]
    mx = max(fin) if fin else None
    print(f"{tag:<24} t={mx if mx else '-':<9} "
          f"{'NORMAL ✅' if 'N O R M A L' in out else 'FAIL ❌'}", flush=True)
