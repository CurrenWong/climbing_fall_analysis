"""快速矩阵：距骨约束 × 多步（落位只用 5 步 ⇒ 单例 ~15 秒）。"""
import os
import re
import subprocess
import sys
from pathlib import Path

R = Path(r"D:\Project\climbing_fall_analysis")
W = R / "temp" / "thums" / "febio_lig"
FEB = r"D:\Program\FEBioStudio\bin\febio4.exe"
PY = str(R / ".venv" / "Scripts" / "python.exe")

CASES = [
    ("transverse + HOLD", dict(FIX_TALUS="transverse", SEAT_HOLD="1")),
    ("transverse 裸", dict(FIX_TALUS="transverse", SEAT_HOLD="0")),
    ("all + HOLD", dict(FIX_TALUS="all", SEAT_HOLD="1")),
    ("none + HOLD", dict(FIX_TALUS="none", SEAT_HOLD="1")),
]
only = sys.argv[1] if len(sys.argv) > 1 else None
print(f"{'case':<24} {'step2→t':>9} {'载荷N':>7}  结果", flush=True)
for tag, ov in CASES:
    if only and only not in tag:
        continue
    e = dict(os.environ)
    e.update(DRIVER="force", FORCE_N="600", TRUNC_MM="60", AXIS_SRC="full",
             ALIGN_AXIS="1", CART_DIR=str(R / "temp" / "thums" / "cart_m1"),
             CART_NEW="1", TIE="1", DET_MIN="0.0", TALUS_RIGID="1",
             TAL_ROTLOCK="0", MULTISTEP="1", SEAT_MM="0.10",
             TIME_STEPS="20", TIME_STEPS_SEAT="5")
    e.update(ov)
    p = subprocess.run([PY, "-u",
                        str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
                       capture_output=True, text=True, env=e, cwd=str(R),
                       timeout=600)
    if p.returncode:
        print(f"{tag:<24}   build 失败 {p.returncode}", flush=True)
        continue
    r = subprocess.run([FEB, "-i", "ankle_lig_WK.feb", "-o", "ms.xplt"],
                       capture_output=True, text=True, cwd=str(W), timeout=900)
    out = (r.stdout or "") + (r.stderr or "")
    fin = [float(x) for x in re.findall(r"converged at time : ([\d.]+)", out)]
    s2 = [t for t in fin if t > 1.0]
    mx = max(s2) if s2 else None
    print(f"{tag:<24} {mx if mx else '-':>9} "
          f"{round(600*(mx-1),1) if mx else 0:>7.1f}  "
          f"{'NORMAL ✅' if 'N O R M A L' in out else 'FAIL ❌'}", flush=True)
