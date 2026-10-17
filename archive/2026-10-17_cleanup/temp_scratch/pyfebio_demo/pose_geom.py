"""验证几何旋转式姿势（POSE_IE=内外翻 / POSE_DF=屈伸）在修复 _outer 后是否可用。"""
import os
import re
import subprocess
from pathlib import Path

R = Path(r"D:\Project\climbing_fall_analysis")
FEB = r"D:\Program\FEBioStudio\bin\febio4.exe"
PY = str(R / ".venv" / "Scripts" / "python.exe")
ACC = str(R / "scripts" / "ankle_fe" / "acceptance.py")
COMMON = ["DRIVER=dispx", "LOAD_SIGN=+1", "FIX_TALUS=transverse", "TRUNC_MM=60",
          "AXIS_SRC=full", "ALIGN_AXIS=1", "CART_NEW=1", "TIE=1", "DET_MIN=0.0",
          "TALUS_RIGID=1", "TAL_ROTLOCK=0", "TIME_STEPS=20", "LAUGON=PENALTY",
          "CART_DIR=" + str(R / "temp" / "thums" / "cart_m1")]
CASES = [
    ("无姿势 DISP=-1.0", ["DISP_MM=-1.0", "POSE_IE=0", "POSE_DF=0"]),
    ("内翻 +10°，DISP=-0.25", ["DISP_MM=-0.25", "POSE_IE=+10"]),
    ("外翻 -10°，DISP=-0.25", ["DISP_MM=-0.25", "POSE_IE=-10"]),
    ("内翻 +10°，DISP=-1.0", ["DISP_MM=-1.0", "POSE_IE=+10"]),
    ("外翻 -10°，DISP=-1.0", ["DISP_MM=-1.0", "POSE_IE=-10"]),
    ("背屈 +10°，DISP=-0.25", ["DISP_MM=-0.25", "POSE_DF=+10"]),
]
print(f"{'case':<24} {'rc':>3} {'NORMAL':>7} {'Fz':>9} {'刚体°':>8} {'内外翻°':>9}",
      flush=True)
for tag, extra in CASES:
    p = subprocess.run([PY, "-u", ACC] + COMMON + extra,
                       capture_output=True, text=True, cwd=str(R), timeout=1800)
    t = (p.stdout or "") + (p.stderr or "")
    ok = "NORMAL TERMINATION" in t
    fz = re.search(r"两顶面合力 = \[\s*[-\d.e+]+\s+[-\d.e+]+\s+([-\d.e+]+)", t)
    rb = re.search(r"全体刚体拟合转动 = ([\d.eE+-]+)", t)
    inv = re.search(r"距骨低端刚体转动 = ([\d.eE+-]+)", t)
    print(f"{tag:<24} {p.returncode:>3} {str(ok):>7} "
          f"{float(fz.group(1)) if fz else 0:>9.2f} "
          f"{float(rb.group(1)) if rb else 0:>8.4f} "
          f"{float(inv.group(1)) if inv else 0:>9.4f}", flush=True)
