"""多步 × 求解器策略（全牛顿 vs BFGS；line search）。"""
import io
import os
import re
import subprocess
from pathlib import Path

R = Path(r"D:\Project\climbing_fall_analysis")
W = R / "temp" / "thums" / "febio_lig"
FEB = r"D:\Program\FEBioStudio\bin\febio4.exe"
CASES = [("基线 BFGS max_ups=10", {}),
         ("全牛顿 max_ups=0", dict(MAX_UPS="0")),
         ("max_ups=5", dict(MAX_UPS="5")),
         ("dtmax 缩小", dict(MAX_DT="0.01"))]
for tag, ov in CASES:
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
    fail = re.findall(r"failed to converge at time : ([\d.]+)", out)
    s2 = [t for t in fin if t > 1.0]
    print(f"{tag:<22} step2 收敛到 t={max(s2) if s2 else '-':<10} "
          f"(载荷 {round(600*(max(s2)-1),1) if s2 else 0:>6} N)  "
          f"失败 t={fail[-1] if fail else '-':<10}  "
          f"{'NORMAL ✅' if 'N O R M A L' in out else 'FAIL ❌'}")
