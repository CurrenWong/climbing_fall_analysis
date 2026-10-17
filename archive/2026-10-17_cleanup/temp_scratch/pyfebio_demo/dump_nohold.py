"""抓「不HOLD + AUGLAG」的失败上下文（是否仍翻负）。"""
import os
import re
import subprocess
from pathlib import Path

R = Path(r"D:\Project\climbing_fall_analysis")
W = R / "temp" / "thums" / "febio_lig"
FEB = r"D:\Program\FEBioStudio\bin\febio4.exe"
PY = str(R / ".venv" / "Scripts" / "python.exe")
e = dict(os.environ)
e.update(DRIVER="force", FORCE_N="600", TRUNC_MM="60", AXIS_SRC="full",
         ALIGN_AXIS="1", CART_DIR=str(R / "temp" / "thums" / "cart_m1"),
         CART_NEW="1", TIE="1", DET_MIN="0.0", TALUS_RIGID="1", TAL_ROTLOCK="0",
         MULTISTEP="1", SEAT_MM="0.10", SEAT_HOLD="0", FIX_TALUS="transverse",
         LAUGON="AUGLAG", TIME_STEPS="20", TIME_STEPS_SEAT="5")
subprocess.run([PY, "-u", str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
               capture_output=True, text=True, env=e, cwd=str(R), timeout=900)
r = subprocess.run([FEB, "-i", "ankle_lig_WK.feb", "-o", "ms.xplt"],
                   capture_output=True, text=True, cwd=str(W), timeout=1800)
out = (r.stdout or "") + (r.stderr or "")
Path(R / "temp" / "pyfebio_demo" / "_k.txt").write_text(out, encoding="utf-8")
lines = [l for l in out.splitlines() if l.strip()]
print("=== 最后 45 行 ===")
print("\n".join(lines[-45:]))
