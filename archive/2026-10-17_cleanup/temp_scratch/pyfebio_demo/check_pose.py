"""查 POSE_TILT 是否真的写进了 BC。"""
import os
import subprocess

R = r"D:\Project\climbing_fall_analysis"
e = dict(os.environ)
e.update(DRIVER="dispx", DISP_MM="-1.0", TRUNC_MM="60", AXIS_SRC="full",
         ALIGN_AXIS="1", CART_DIR=R + r"\temp\thums\cart_m1", CART_NEW="1",
         TIE="1", DET_MIN="0.0", TALUS_RIGID="1", TAL_ROTLOCK="0",
         FIX_TALUS="transverse", TIME_STEPS="20", POSE_TILT="10")
p = subprocess.run([R + r"\.venv\Scripts\python.exe", "-u",
                    R + r"\scripts\ankle_fe\ankle_lig_feb.py"],
                   capture_output=True, text=True, env=e, cwd=R, timeout=900)
t = (p.stdout or "") + (p.stderr or "")
for l in t.splitlines():
    if "pose" in l.lower() or "倾斜" in l:
        print("LOG:", l)
print("rc", p.returncode)
