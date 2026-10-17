"""夹出 600 N 对应的位移（力的目标值可达性检查）。"""
import os
import re
import subprocess

ROOT = r"D:\Project\climbing_fall_analysis"
PY = rf"{ROOT}\.venv\Scripts\python.exe"
ACC = rf"{ROOT}\scripts\ankle_fe\acceptance.py"
BASE = dict(DRIVER="dispx", FIX_TALUS="transverse", TRUNC_MM="60", AXIS_SRC="full",
            ALIGN_AXIS="1", CART_DIR=rf"{ROOT}\temp\thums\cart_m1", CART_NEW="1",
            TIE="1", DET_MIN="0.0")
print(f"{'DISP(mm)':>9} {'状态':>7} {'Fz(N)':>10} {'k=F/δ':>9}  {'|u|max':>9}")
print("-" * 52)
for disp in ("1.50", "1.60", "1.70", "1.80"):
    e = dict(os.environ)
    e.update(BASE)
    e["DISP_MM"] = disp
    r = subprocess.run([PY, "-u", ACC] + [f"{k}={v}" for k, v in
                                          list(BASE.items()) + [("DISP_MM", disp)]],
                       capture_output=True, text=True, env=e, cwd=ROOT, timeout=1800)
    t = r.stdout + r.stderr
    norm = "NORMAL" if "NORMAL TERMINATION" in t else ("ERROR" if "ERROR" in t else "?")
    fz = re.search(r"两顶面合力 = \[\s*[-\d.e+]+\s+[-\d.e+]+\s+([-\d.e+]+)\]", t)
    um = re.search(r"\|u\|max = ([\d.]+)", t)
    v = float(fz.group(1)) if fz else float("nan")
    print(f"{disp:>9} {norm:>7} {v:>10.2f} {v/float(disp):>9.1f}  "
          f"{(um.group(1) if um else '-'):>9}")
