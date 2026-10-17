"""路线 A 二分：逐步关掉 tie / 软骨新格式，定位 rigid 初始化失败的原因。"""
import os
import re
import subprocess

ROOT = r"D:\Project\climbing_fall_analysis"
PY = rf"{ROOT}\.venv\Scripts\python.exe"
SCRIPT = rf"{ROOT}\scripts\ankle_fe\ankle_lig_feb.py"
BASE = dict(DRIVER="dispx", DISP_MM="0.25", FIX_TALUS="transverse", TRUNC_MM="60",
            AXIS_SRC="full", ALIGN_AXIS="1", TALUS_RIGID="1", DET_MIN="0.0")

CASES = [
    ("① 旧软骨 + 无 tie", dict(TIE="0")),
    ("② m1 新软骨 + 无 tie", dict(CART_DIR=rf"{ROOT}\temp\thums\cart_m1",
                                CART_NEW="1", TIE="0")),
    ("③ m1 新软骨 + tie", dict(CART_DIR=rf"{ROOT}\temp\thums\cart_m1",
                              CART_NEW="1", TIE="1")),
]
for tag, ov in CASES:
    e = dict(os.environ)
    e.update(BASE)
    e.update(ov)
    r = subprocess.run([PY, "-u", SCRIPT], capture_output=True, text=True,
                       env=e, cwd=ROOT, timeout=2400)
    t = r.stdout + r.stderr
    rc = re.findall(r"FEBio rc=(\d+)", t)
    bad = [l.strip() for l in t.splitlines()
           if "onstraint" in l or "NORMAL" in l or "tag " in l
           or "unrecognized" in l or "invalid" in l]
    print(f"\n===== {tag} =====  rc={rc}")
    for l in bad[:6]:
        print("   " + l[:100])
