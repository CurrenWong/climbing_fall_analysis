"""抓多步运行的 FEBio 输出（step1 是否完成、step2 为何失败）。"""
import io
import os
import re
import subprocess

R = r"D:\Project\climbing_fall_analysis"
LOG = os.path.join(R, "temp", "pyfebio_demo", "_ms_febio.txt")
e = dict(os.environ)
e.update(DRIVER="force", FORCE_N="600", FIX_TALUS="transverse", TRUNC_MM="60",
         AXIS_SRC="full", ALIGN_AXIS="1", CART_DIR=R + r"\temp\thums\cart_m1",
         CART_NEW="1", TIE="1", DET_MIN="0.0", TALUS_RIGID="1",
         TAL_ROTLOCK="0", MULTISTEP="1", SEAT_MM="0.5")
r = subprocess.run([R + r"\.venv\Scripts\python.exe", "-u",
                    R + r"\scripts\ankle_fe\ankle_lig_feb.py"],
                   capture_output=True, text=True, env=e, cwd=R, timeout=3000)
t = r.stdout + r.stderr
# 找 step 相关与失败上下文
keys = ("beginning time step", "step", "converged", "failed", "ERROR",
        "negative", "reformation", "Augmentation", "residual")
lines = t.splitlines()
idx = [i for i, l in enumerate(lines) if "failed" in l or "ERROR" in l]
k = idx[0] if idx else len(lines) - 1
print("=== 失败前 30 行 ===")
print("\n".join(lines[max(0, k - 30):k + 3]))
print("\n=== t≈1.0 附近的步进记录 ===")
for i, l in enumerate(lines):
    if re.search(r"beginning time step|converged at time", l):
        print("  " + l.strip()[:80])
