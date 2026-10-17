"""细扫 TIE 罚因子（auto=0），找力是否收敛到平台。"""
import os
import re
import subprocess

ROOT = r"D:\Project\climbing_fall_analysis"
PY = rf"{ROOT}\.venv\Scripts\python.exe"
ACC = rf"{ROOT}\scripts\ankle_fe\acceptance.py"
BASE = dict(DRIVER="dispx", DISP_MM="1.00", FIX_TALUS="transverse", TRUNC_MM="60",
            AXIS_SRC="full", ALIGN_AXIS="1",
            CART_DIR=rf"{ROOT}\temp\thums\cart_m1", CART_NEW="1", TIE="1", DET_MIN="0.0")

print(f"{'TIE_PENALTY (auto=0)':<24} {'状态':>7} {'Fz(N)':>11}")
print("-" * 46)
for pen in ("0.1", "1", "3", "10", "30", "100", "300"):
    ov = {"TIE_PENALTY": pen, "TIE_AUTO_PENALTY": "0"}
    e = dict(os.environ)
    e.update(BASE)
    e.update(ov)
    r = subprocess.run([PY, "-u", ACC] + [f"{k}={v}" for k, v in
                                          list(BASE.items()) + list(ov.items())],
                       capture_output=True, text=True, env=e, cwd=ROOT, timeout=1800)
    t = r.stdout + r.stderr
    norm = "NORMAL" if "NORMAL TERMINATION" in t else ("ERROR" if "ERROR" in t else "?")
    fz = re.search(r"两顶面合力 = \[\s*[-\d.e+]+\s+[-\d.e+]+\s+([-\d.e+]+)\]", t)
    print(f"{pen:<24} {norm:>7} {(float(fz.group(1)) if fz else float('nan')):>11.2f}")
