"""多步落位首次运行。"""
import os
import re
import subprocess

R = r"D:\Project\climbing_fall_analysis"
e = dict(os.environ)
e.update(DRIVER="force", FORCE_N="600", FIX_TALUS="transverse", TRUNC_MM="60",
         AXIS_SRC="full", ALIGN_AXIS="1", CART_DIR=R + r"\temp\thums\cart_m1",
         CART_NEW="1", TIE="1", DET_MIN="0.0", TALUS_RIGID="1",
         TAL_ROTLOCK="0", MULTISTEP="1", SEAT_MM="0.5")
r = subprocess.run([R + r"\.venv\Scripts\python.exe", "-u",
                    R + r"\scripts\ankle_fe\ankle_lig_feb.py"],
                   capture_output=True, text=True, env=e, cwd=R, timeout=3000)
t = r.stdout + r.stderr
print("multistep:", re.findall(r"\[multistep\].*", t))
print("talus-rigid:", re.findall(r"\[talus-rigid\].*", t)[:1])
print("rc:", re.findall(r"FEBio rc=.*", t))
print("err:")
for l in t.splitlines():
    if any(k in l for k in ("tag ", "invalid", "unrecognized", "ERROR",
                            "onverge", "negative")):
        print("   " + l.strip()[:105])
print("ok:", [l.strip() for l in t.splitlines()
              if "converged" in l or "NORMAL" in l][:4])
