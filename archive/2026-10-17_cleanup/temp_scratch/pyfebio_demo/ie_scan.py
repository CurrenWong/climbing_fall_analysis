"""看几何姿势的落位收敛 + 扫位移找稳定区。"""
import os
import re
import subprocess
from pathlib import Path

R = Path(r"D:\Project\climbing_fall_analysis")
PY = str(R / ".venv" / "Scripts" / "python.exe")
ACC = str(R / "scripts" / "ankle_fe" / "acceptance.py")
COMMON = ["DRIVER=dispx", "LOAD_SIGN=+1", "FIX_TALUS=transverse", "TRUNC_MM=60",
          "AXIS_SRC=full", "ALIGN_AXIS=1", "CART_NEW=1", "TIE=1", "DET_MIN=0.0",
          "TALUS_RIGID=1", "TAL_ROTLOCK=0", "TIME_STEPS=20", "LAUGON=PENALTY",
          "CART_DIR=" + str(R / "temp" / "thums" / "cart_m1")]


def go(tag, extra, show=False):
    p = subprocess.run([PY, "-u", ACC] + COMMON + extra,
                       capture_output=True, text=True, cwd=str(R), timeout=1800)
    t = (p.stdout or "") + (p.stderr or "")
    if show:
        for l in t.splitlines():
            if "[pose]" in l or "负 Jacobian" in l or "求解" in l:
                print("   ", l.strip())
    ok = "NORMAL TERMINATION" in t
    fz = re.search(r"两顶面合力 = \[\s*[-\d.e+]+\s+[-\d.e+]+\s+([-\d.e+]+)", t)
    return ok, (float(fz.group(1)) if fz else 0.0)


print("--- 落位明细（POSE_IE=+10, DISP=-0.25）---")
go("x", ["DISP_MM=-0.25", "POSE_IE=+10"], show=True)

print("\n--- 扫 DISP_MM（POSE_IE=+10 内翻 / −10 外翻）---")
for ie in ("+10", "-10"):
    for d in ("-0.25", "-0.5", "-0.8", "-1.0", "-1.2"):
        ok, fz = go("x", [f"DISP_MM={d}", f"POSE_IE={ie}"])
        print(f"  POSE_IE={ie:>4} DISP={d:>6}  {'NORMAL ✅' if ok else 'FAIL ❌'}"
              f"  Fz={fz:>9.2f}", flush=True)
