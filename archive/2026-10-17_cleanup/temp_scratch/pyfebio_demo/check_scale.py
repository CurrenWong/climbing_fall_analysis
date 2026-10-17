"""确认 <scale lc> 已写入，并看 step2 失败的当前模式。"""
import io
import re
import subprocess
from pathlib import Path

R = Path(r"D:\Project\climbing_fall_analysis")
W = R / "temp" / "thums" / "febio_lig"
FEB = r"D:\Program\FEBioStudio\bin\febio4.exe"
import os
e = dict(os.environ)
e.update(DRIVER="force", FORCE_N="600", FIX_TALUS="transverse", TRUNC_MM="60",
         AXIS_SRC="full", ALIGN_AXIS="1",
         CART_DIR=str(R / "temp" / "thums" / "cart_m1"), CART_NEW="1", TIE="1",
         DET_MIN="0.0", TALUS_RIGID="1", TAL_ROTLOCK="0", MULTISTEP="1",
         SEAT_MM="0.10", TIME_STEPS="40", TIME_STEPS_SEAT="20")
subprocess.run([str(R / ".venv" / "Scripts" / "python.exe"), "-u",
                str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
               capture_output=True, text=True, env=e, cwd=str(R), timeout=3000)
f = (W / "ankle_lig_WK.feb").read_text(encoding="latin-1", errors="replace")
k = f.find('<step id="2">')
print("=== step2 ===")
print(f[k:k + 1400])
