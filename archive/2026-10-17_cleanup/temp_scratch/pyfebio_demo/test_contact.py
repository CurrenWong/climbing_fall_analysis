"""多步 × 接触算法（LAUGON / NODE_RELOC）。"""
import os
import re
import subprocess
from pathlib import Path

R = Path(r"D:\Project\climbing_fall_analysis")
FEB = r"D:\Program\FEBioStudio\bin\febio4.exe"
W = R / "temp" / "thums" / "febio_lig"
CASES = [
    ("PENALTY  seat=0.5", dict(LAUGON="PENALTY", SEAT_MM="0.5")),
    ("AUGLAG   seat=0.5", dict(LAUGON="AUGLAG", SEAT_MM="0.5")),
    ("AUGLAG   seat=0.9", dict(LAUGON="AUGLAG", SEAT_MM="0.9")),
    ("AUGLAG+reloc seat=0.5", dict(LAUGON="AUGLAG", SEAT_MM="0.5",
                                   NODE_RELOC="1")),
    ("PENALTY+reloc seat=0.5", dict(LAUGON="PENALTY", SEAT_MM="0.5",
                                    NODE_RELOC="1")),
]
for tag, ov in CASES:
    e = dict(os.environ)
    e.update(DRIVER="force", FORCE_N="600", FIX_TALUS="transverse", TRUNC_MM="60",
             AXIS_SRC="full", ALIGN_AXIS="1",
             CART_DIR=str(R / "temp" / "thums" / "cart_m1"), CART_NEW="1",
             TIE="1", DET_MIN="0.0", TALUS_RIGID="1", TAL_ROTLOCK="0",
             MULTISTEP="1", TIME_STEPS="20", TIME_STEPS_SEAT="20")
    e.update(ov)
    subprocess.run([str(R / ".venv" / "Scripts" / "python.exe"), "-u",
                    str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
                   capture_output=True, text=True, env=e, cwd=str(R), timeout=3000)
    r = subprocess.run([FEB, "-i", "ankle_lig_WK.feb", "-o", "ms.xplt"],
                       capture_output=True, text=True, cwd=str(W), timeout=3000)
    out = (r.stdout or "") + (r.stderr or "")
    neg = re.findall(r"(\d+) negative jacobian", out)
    fail = re.findall(r"failed to converge at time : ([\d.]+)", out)
    fin = re.findall(r"converged at time : ([\d.]+)", out)
    norm = "N O R M A L" in out
    print(f"{tag:<24} 收敛末 t={fin[-1] if fin else '-':<10} "
          f"失败 t={fail[-1] if fail else '-':<10} negJ={neg[:1]}  "
          f"{'NORMAL ✅' if norm else 'FAIL ❌'}")
