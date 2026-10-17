"""干净软骨网格(cart_patch) + TIE 的力-位移曲线。

对照目标：旧格式软骨在**位移控制**下 ~0.3mm 就翻负（≈235 N 天花板）。
本脚本逐点跑，读回 两顶面合力 与 |u|max。
"""
import os
import re
import subprocess

ROOT = r"D:\Project\climbing_fall_analysis"
PY = rf"{ROOT}\.venv\Scripts\python.exe"
ACC = rf"{ROOT}\scripts\ankle_fe\acceptance.py"

BASE = dict(DRIVER="dispx", FIX_TALUS="transverse", TRUNC_MM="60", AXIS_SRC="full",
            ALIGN_AXIS="1", CART_DIR=rf"{ROOT}\temp\thums\cart_m1", CART_NEW="1",
            TIE="1", DET_MIN="0.0")

print(f"{'DISP(mm)':>9} {'状态':>6} {'合力Fz(N)':>12} {'|u|max(mm)':>11}  备注")
print("-" * 62)
for disp in ("0.10", "0.25", "0.50", "0.75", "1.00", "1.50", "2.00"):
    e = dict(os.environ)
    e.update(BASE)
    e["DISP_MM"] = disp
    r = subprocess.run([PY, "-u", ACC] + [f"{k}={v}" for k, v in
                                          list(BASE.items()) + [("DISP_MM", disp)]],
                       capture_output=True, text=True, env=e, cwd=ROOT, timeout=1800)
    t = r.stdout + r.stderr
    norm = "NORMAL" if "NORMAL TERMINATION" in t else (
        "ERROR" if "ERROR" in t else "?")
    fj = re.search(r"两顶面合力 = \[[^]]*\]\s*N\s*\|\s*F\|=([\d.]+)", t)
    um = re.search(r"\|u\|max = ([\d.]+)", t)
    fz = re.search(r"两顶面合力 = \[\s*([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)\]", t)
    fzs = fz.group(3) if fz else "-"
    print(f"{disp:>9} {norm:>6} {str(fzs):>12} {um.group(1) if um else '-':>11}")
