"""零回归（正确传参）：MULTISTEP 默认关 ⇒ 基线数值应逐字不变。"""
import os
import re
import subprocess

R = r"D:\Project\climbing_fall_analysis"
PY = R + r"\.venv\Scripts\python.exe"
COMMON = ["DET_MIN=0.0", "TRUNC_MM=60", "AXIS_SRC=full", "ALIGN_AXIS=1",
          "FIX_TALUS=transverse"]
CASES = [
    ("默认路径(旧格式)", ["DRIVER=disp", "DISP_MM=0.1472"]),
    ("TIE=1+CART_NEW m=1", ["DRIVER=disp", "DISP_MM=1.6", "TIE=1", "CART_NEW=1",
                            "CART_DIR=" + R + r"\temp\thums\cart_m1"]),
    ("路线A rigid 位移控制", ["TALUS_RIGID=1", "DRIVER=dispx", "DISP_MM=0.25"]),
    ("路线A 转轴解锁", ["TALUS_RIGID=1", "DRIVER=dispx", "DISP_MM=0.25",
                    "TAL_ROTLOCK=0"]),
]
for tag, args in CASES:
    r = subprocess.run([PY, "-u", R + r"\scripts\ankle_fe\acceptance.py"]
                       + COMMON + args,
                       capture_output=True, text=True, cwd=R, timeout=3000)
    t = r.stdout + r.stderr
    nm = ("NORMAL" if "NORMAL TERMINATION" in t
          else ("ERROR" if "E R R O R" in t else "?"))
    fz = re.search(r"两顶面合力 = \[\s*[-\d.e+]+\s+[-\d.e+]+\s+([-\d.e+]+)", t)
    print(f"{tag:<22} {nm:>6}  Fz={float(fz.group(1)) if fz else float('nan'):>10.4f} N")
