"""聚焦批量：围绕"单步 + 压缩 + PENALTY ts=20"（已到 337 N）找能通关的组合。"""
import os
import re
import subprocess
from pathlib import Path

R = Path(r"D:\Project\climbing_fall_analysis")
W = R / "temp" / "thums" / "febio_lig"
FEB = r"D:\Program\FEBioStudio\bin\febio4.exe"
PY = str(R / ".venv" / "Scripts" / "python.exe")
CASES = [
    ("基线 ts=20 P", {}),
    ("ts=20 P NODE_RELOC", dict(NODE_RELOC="1")),
    ("ts=20 P 变形距骨", dict(TALUS_RIGID="0")),
    ("ts=20 P 转轴锁", dict(TAL_ROTLOCK="1")),
    ("ts=20 P penalty=10", dict(PENALTY="10")),
    ("ts=20 P penalty=0.1", dict(PENALTY="0.1")),
]
for tag, ov in CASES:
    e = dict(os.environ)
    e.update(DRIVER="force", FORCE_N="600", TRUNC_MM="60", AXIS_SRC="full",
             ALIGN_AXIS="1", CART_DIR=str(R / "temp" / "thums" / "cart_m1"),
             CART_NEW="1", TIE="1", DET_MIN="0.0", TALUS_RIGID="1",
             TAL_ROTLOCK="0", FIX_TALUS="transverse", TIME_STEPS="20",
             LOAD_SIGN="+1", LAUGON="PENALTY")
    e.update(ov)
    subprocess.run([PY, "-u",
                    str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
                   capture_output=True, text=True, env=e, cwd=str(R), timeout=900)
    r = subprocess.run([FEB, "-i", "ankle_lig_WK.feb", "-o", "s.xplt"],
                       capture_output=True, text=True, cwd=str(W), timeout=1800)
    out = (r.stdout or "") + (r.stderr or "")
    fin = [float(x) for x in re.findall(r"converged at time : ([\d.]+)", out)]
    mx = max(fin) if fin else None
    print(f"{tag:<22} t={mx if mx else '-':<9} "
          f"载荷≈{round(600*mx,0) if mx else 0:>5.0f} N  "
          f"{'NORMAL ✅' if 'N O R M A L' in out else 'FAIL ❌'}", flush=True)
