"""看 POSE_IE=+10 构建/求解的失败细节。"""
import subprocess

R = r"D:\Project\climbing_fall_analysis"
A = [R + r"\.venv\Scripts\python.exe", "-u",
     R + r"\scripts\ankle_fe\acceptance.py", "DRIVER=dispx", "DISP_MM=-0.25",
     "LOAD_SIGN=+1", "FIX_TALUS=transverse", "TRUNC_MM=60", "AXIS_SRC=full",
     "ALIGN_AXIS=1", "CART_NEW=1", "TIE=1", "DET_MIN=0.0", "TALUS_RIGID=1",
     "TAL_ROTLOCK=0", "TIME_STEPS=20", "LAUGON=PENALTY", "POSE_IE=+10",
     "CART_DIR=" + R + r"\temp\thums\cart_m1"]
p = subprocess.run(A, capture_output=True, text=True, cwd=R, timeout=1800)
t = (p.stdout or "") + (p.stderr or "")
print(t[:3200])
