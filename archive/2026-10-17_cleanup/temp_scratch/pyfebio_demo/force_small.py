"""小微载力控制测试：验证「卡在关节初始张开」的假设。"""
import os
import re
import subprocess

R = r"D:\Project\climbing_fall_analysis"
PY = R + r"\.venv\Scripts\python.exe"
ACC = R + r"\scripts\ankle_fe\acceptance.py"
BASE = dict(DRIVER="force", FIX_TALUS="transverse", TRUNC_MM="60", AXIS_SRC="full",
            ALIGN_AXIS="1", CART_DIR=R + r"\temp\thums\cart_m1", CART_NEW="1",
            TIE="1", DET_MIN="0.0", TALUS_RIGID="1", TAL_ROTLOCK="0")
for fn in ("1", "5", "20", "60"):
    d = dict(BASE)
    d["FORCE_N"] = fn
    e = dict(os.environ)
    e.update(d)
    r = subprocess.run([PY, "-u", ACC] + [f"{k}={v}" for k, v in d.items()],
                       capture_output=True, text=True, env=e, cwd=R, timeout=2400)
    t = r.stdout + r.stderr
    nm = ("NORMAL" if "NORMAL TERMINATION" in t else
          ("ERROR" if "E R R O R" in t else "?"))
    nj = re.findall(r"(\d+) negative jacobian", t)
    um = re.search(r"\|u\|max = ([\d.]+)", t)
    print(f"  FORCE_N={fn:>3}  {nm:>6}  negJ={nj[:1]}  |u|max={um.group(1) if um else '-'}")
