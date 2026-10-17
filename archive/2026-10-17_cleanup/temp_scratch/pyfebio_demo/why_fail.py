"""看 build 失败原因。"""
import os
import subprocess

R = r"D:\Project\climbing_fall_analysis"
e = dict(os.environ)
e.update(DRIVER="force", FORCE_N="600", TRUNC_MM="60", AXIS_SRC="full",
         ALIGN_AXIS="1", CART_DIR=R + r"\temp\thums\cart_m1", CART_NEW="1",
         TIE="1", DET_MIN="0.0", TALUS_RIGID="1", TAL_ROTLOCK="0",
         FIX_TALUS="transverse", TIME_STEPS="20", LAUGON="AUGLAG",
         LOAD_SIGN="+1")
p = subprocess.run([R + r"\.venv\Scripts\python.exe", "-u",
                    R + r"\scripts\ankle_fe\ankle_lig_feb.py"],
                   capture_output=True, text=True, env=e, cwd=R, timeout=900)
t = (p.stdout or "") + (p.stderr or "")
print("rc", p.returncode)
print(t[-2500:])
