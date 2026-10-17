"""A4：用更细的软骨网格（cart_m3，m=3）测 A1+A2 的量程。

对照基线：cart_m1（m=1，182 个胫距距骨单元）
测：cart_m3 在 +1.6 / +1.4 / +1.2 / +1.0 下是否收敛
"""
import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WORK = ROOT / "temp/thums/febio_lig"
OUT = ROOT / "temp/pyfebio_demo/a4mesh"
OUT.mkdir(parents=True, exist_ok=True)
PY = ROOT / ".venv/Scripts/python.exe"
SCRIPT = ROOT / "scripts/ankle_fe/ankle_lig_feb.py"
BASE = dict(DRIVER="dispx", FIX_TALUS="transverse", DET_MIN="0.0",
            TRUNC_MM="60", AXIS_SRC="full", ALIGN_AXIS="1",
            CART_NEW="1", TIE="1", POSE_IE="0", LOG_DATA="1",
            SUBTALAR_LAUGON="PENALTY", SUBTALAR="1",
            TALUS_SUPPORT="spring", TALUS_K="500", MAX_UPS="0")


def sum_top(p):
    if not p.exists():
        return 0.0
    cur_t, buf, last = None, [], {}
    for l in p.read_text(encoding="utf-8", errors="replace").split("\n"):
        if l.startswith("*Time"):
            if buf:
                last[cur_t] = np.array(buf, float)[:, 1:4].sum(0)
            cur_t = float(l.split("=")[1]); buf = []
        elif l.strip() and not l.startswith("*") and not l.startswith("*Step"):
            try:
                buf.append([float(x) for x in l.split(",")])
            except ValueError:
                pass
    if buf:
        last[cur_t] = np.array(buf, float)[:, 1:4].sum(0)
    if not last:
        return 0.0
    return float(last[max(last)][2])


def run(tag, cart, extra):
    env = os.environ.copy()
    env.update(BASE); env.update({"CART_DIR": str(ROOT / "temp/thums" / cart)})
    env.update(extra); env["PYTHONIOENCODING"] = "utf-8"
    try:
        r = subprocess.run([str(PY), "-u", str(SCRIPT)], cwd=str(ROOT), env=env,
                           capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=2400)
    except subprocess.TimeoutExpired:
        print(f"[{tag}] 超时"); return False, None
    out = (r.stdout or "") + (r.stderr or "")
    (OUT / f"{tag}.stdout.txt").write_text(out, encoding="utf-8")
    lg = WORK / "log_WK.txt"
    normal = lg.exists() and "N O R M A L   T E R M I N A T I O N" in \
        lg.read_text(encoding="utf-8", errors="replace")
    for src, dst in ((WORK / "febio_WK.log", f"{tag}.febio.log"),
                     (WORK / "ankle_lig_WK.feb", f"{tag}.feb")):
        if src.exists():
            shutil.copy(src, OUT / dst)
    F = (sum_top(WORK / "reac_tibia_WK.txt")
         + sum_top(WORK / "reac_fibula_WK.txt")) if normal else None
    print(f"[{tag}] cart={cart} NORMAL={str(normal):<5} F={F}")
    return normal, F


print("=== cart_m3（m=3）===")
for d in ("+1.6", "+1.4", "+1.2", "+1.0"):
    run(f"m3_{d}", "cart_m3", {"DISP_MM": d})

print("\n=== 对照：cart_m1（m=1）在 +1.2 / +1.0 ===")
for d in ("+1.2", "+1.0"):
    run(f"m1_{d}", "cart_m1", {"DISP_MM": d})
