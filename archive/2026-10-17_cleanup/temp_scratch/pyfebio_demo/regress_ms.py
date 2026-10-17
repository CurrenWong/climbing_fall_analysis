"""零回归：MULTISTEP 默认关，基线数值应逐字不变。"""
import os
import re
import subprocess

R = r"D:\Project\climbing_fall_analysis"
PY = R + r"\.venv\Scripts\python.exe"
CASES = [
    ("默认路径(旧格式)", dict(DISP_MM="0.1472")),
    ("TIE=1+CART_NEW m=1", dict(TIE="1", CART_NEW="1",
                              CART_DIR=R + r"\temp\thums\cart_m1", DISP_MM="1.6")),
    ("路线A rigid 位移控制", dict(TALUS_RIGID="1", DRIVER="dispx", DISP_MM="0.25")),
]
for tag, ov in CASES:
    e = dict(os.environ)
    e.update(FIX_TALUS="transverse", TRUNC_MM="60", AXIS_SRC="full", ALIGN_AXIS="1",
             DET_MIN="0.0")
    e.update(ov)
    r = subprocess.run([PY, "-u", R + r"\scripts\ankle_fe\acceptance.py"] +
                       [f"{k}={v}" for k, v in e.items() if k in ov or k in
                        ("DRIVER", "DISP_MM", "FIX_TALUS", "TRUNC_MM", "AXIS_SRC",
                         "ALIGN_AXIS", "DET_MIN", "TIE", "CART_NEW", "CART_DIR",
                         "TALUS_RIGID") and k in e],
                       capture_output=True, text=True, env=e, cwd=R, timeout=3000)
    t = r.stdout + r.stderr
    nm = ("NORMAL" if "NORMAL TERMINATION" in t
          else ("ERROR" if "E R R O R" in t else "?"))
    fz = re.search(r"两顶面合力 = \[\s*[-\d.e+]+\s+[-\d.e+]+\s+([-\d.e+]+)", t)
    print(f"{tag:<22} {nm:>6}  Fz={float(fz.group(1)) if fz else float('nan'):>10.4f} N")
