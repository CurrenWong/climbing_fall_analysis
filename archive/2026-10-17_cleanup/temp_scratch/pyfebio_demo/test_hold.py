"""SEAT_HOLD 对照：不撤约束时载荷能否叠加（检验"回弹"机制）。"""
import os
import re
import subprocess
from pathlib import Path

R = Path(r"D:\Project\climbing_fall_analysis")
W = R / "temp" / "thums" / "febio_lig"
FEB = r"D:\Program\FEBioStudio\bin\febio4.exe"
for tag, ov in (("撤约束 (基线)", dict(SEAT_HOLD="0")),
                ("保持约束 SEAT_HOLD=1", dict(SEAT_HOLD="1"))):
    e = dict(os.environ)
    e.update(DRIVER="force", FORCE_N="600", FIX_TALUS="transverse", TRUNC_MM="60",
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
    s2 = [t for t in fin if t > 1.0]
    neg = re.findall(r"(\d+) negative jacobian", out)
    print(f"{tag:<24} step2 到 t={max(s2) if s2 else '-':<10} "
          f"(载荷≈{round(600*(max(s2)-1),1) if s2 else 0:>6} N) negJ={neg[:1]} "
          f"{'NORMAL ✅' if 'N O R M A L' in out else 'FAIL ❌'}")
