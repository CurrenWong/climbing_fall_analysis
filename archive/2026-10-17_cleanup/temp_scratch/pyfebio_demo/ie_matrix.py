"""真实内外翻姿势矩阵（几何旋转式）—— 判据 C2（非零内外翻角）+ 载荷 + 力臂效应。"""
import os
import re
import subprocess
from pathlib import Path

R = Path(r"D:\Project\climbing_fall_analysis")
PY = str(R / ".venv" / "Scripts" / "python.exe")
ACC = str(R / "scripts" / "ankle_fe" / "acceptance.py")
COMMON = ["DRIVER=dispx", "DISP_MM=-1.0", "LOAD_SIGN=+1",
          "FIX_TALUS=transverse", "TRUNC_MM=60", "AXIS_SRC=full",
          "ALIGN_AXIS=1", "CART_NEW=1", "TIE=1", "DET_MIN=0.0",
          "TALUS_RIGID=0", "TAL_ROTLOCK=0", "TIME_STEPS=20",
          "LAUGON=PENALTY",
          "CART_DIR=" + str(R / "temp" / "thums" / "cart_m1")]
print(f"{'POSE_IE':>8} {'PASS(C1)':>9} {'Fz':>9} {'Fx':>8} {'|F|':>8} "
      f"{'刚体°':>8} {'距骨转动°':>10}", flush=True)
rows = []
for ie in ("-5", "-3", "-1", "0", "+1", "+3", "+5"):
    p = subprocess.run([PY, "-u", ACC] + COMMON + [f"POSE_IE={ie}"],
                       capture_output=True, text=True, cwd=str(R), timeout=1800)
    t = (p.stdout or "") + (p.stderr or "")
    ok = "NORMAL TERMINATION" in t
    v = re.search(r"两顶面合力 = \[\s*([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)", t)
    rb = re.search(r"全体刚体拟合转动 = ([\d.eE+-]+)", t)
    tl = re.search(r"距骨低端刚体转动 = ([\d.eE+-]+)", t)
    fx, fy, fz = (float(v.group(1)), float(v.group(2)), float(v.group(3))) \
        if v else (0, 0, 0)
    rows.append((ie, ok, fx, fz, rb.group(1) if rb else "?", tl.group(1) if tl else "?"))
    print(f"{ie:>8} {str(ok):>9} {fz:>9.2f} {fx:>8.2f} "
          f"{(fx*fx+fz*fz)**0.5:>8.2f} "
          f"{float(rb.group(1)) if rb else 0:>8.4f} "
          f"{float(tl.group(1)) if tl else 0:>10.4f}", flush=True)

print("\n=== C2（非零内外翻角）===")
for ie, ok, fx, fz, rb, tl in rows:
    verdict = "PASS" if (ok and float(tl) > 0.05) else "FAIL"
    print(f"  POSE_IE={ie:>3}  NORMAL={ok}  距骨刚体转动={tl}°  ⇒ C2 {verdict}")
