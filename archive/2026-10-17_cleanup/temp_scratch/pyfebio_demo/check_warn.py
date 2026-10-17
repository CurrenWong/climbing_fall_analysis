"""确认告警在默认（拉伸）方向下会打印，且压缩方向不打印。"""
import os
import subprocess

R = r"D:\Project\climbing_fall_analysis"
PY = R + r"\.venv\Scripts\python.exe"
SCR = R + r"\scripts\ankle_fe\ankle_lig_feb.py"
BASE = dict(DRIVER="dispx", DISP_MM="-1.0", TRUNC_MM="60", AXIS_SRC="full",
            ALIGN_AXIS="1", CART_DIR=R + r"\temp\thums\cart_m1", CART_NEW="1",
            TIE="1", DET_MIN="0.0", TALUS_RIGID="1", TAL_ROTLOCK="0",
            FIX_TALUS="transverse", TIME_STEPS="20")
for tag, ov in (("默认(拉伸)", {}), ("压缩 LOAD_SIGN=+1", dict(LOAD_SIGN="+1"))):
    e = dict(os.environ)
    e.update(BASE)
    e.update(ov)
    p = subprocess.run([PY, "-u", SCR], capture_output=True, text=True, env=e,
                       cwd=R, timeout=900)
    print(f"--- {tag} (rc={p.returncode}) ---")
    for l in ((p.stdout or "") + (p.stderr or "")).splitlines():
        if "[load]" in l:
            print("  " + l)
