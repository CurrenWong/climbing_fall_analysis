"""跑 m=1 并完整打印结果。"""
import os
import re
import subprocess

ROOT = r"D:\Project\climbing_fall_analysis"
for m in ("1", "3"):
    e = dict(os.environ)
    e.update(dict(DRIVER="dispx", DISP_MM="0.25", FIX_TALUS="transverse", TRUNC_MM="60",
                  AXIS_SRC="full", ALIGN_AXIS="1",
                  CART_DIR=rf"{ROOT}\temp\thums\cart_m{m}", CART_NEW="1", TIE="1",
                  DET_MIN="0.0"))
    r = subprocess.run([rf"{ROOT}\.venv\Scripts\python.exe", "-u",
                        rf"{ROOT}\scripts\ankle_fe\ankle_lig_feb.py"],
                       capture_output=True, text=True, env=e, cwd=ROOT)
    t = r.stdout + r.stderr
    print(f"\n{'='*70}\n### m={m}   rc={r.returncode}\n{'='*70}")
    keys = ("[cart]", "[tie]", "接触面", "翻转修正", "negative jacobian",
            "Negative jacobian", "converge", "NORMAL", "E R R O R", "写出")
    for ln in t.splitlines():
        if any(k in ln for k in keys):
            print("  " + ln.strip())
