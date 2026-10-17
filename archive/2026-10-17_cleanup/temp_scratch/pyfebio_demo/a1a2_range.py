"""A1+A2 载荷量程测试：E2p 配置在 DISP=+1.6/+1.0/+0.5/+0.25 下是否收敛。
目的：区分「载荷过高」vs「配置本身坏了」，并给出可用的位移量程。
对照：同一载荷下基线（SUBTALAR=0, fix）都应收敛。
"""
import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WORK = ROOT / "temp/thums/febio_lig"
OUT = ROOT / "temp/pyfebio_demo/a1a2"
PY = ROOT / ".venv/Scripts/python.exe"
SCRIPT = ROOT / "scripts/ankle_fe/ankle_lig_feb.py"
CART = str(ROOT / "temp/thums/cart_m1")
BASE = dict(DRIVER="dispx", FIX_TALUS="transverse", DET_MIN="0.0", TRUNC_MM="60",
            AXIS_SRC="full", ALIGN_AXIS="1", CART_NEW="1", TIE="1",
            CART_DIR=CART, POSE_IE="0", LOG_DATA="1", SUBTALAR_LAUGON="PENALTY")


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


def run(tag, extra):
    env = os.environ.copy()
    env.update(BASE); env.update(extra); env["PYTHONIOENCODING"] = "utf-8"
    try:
        r = subprocess.run([str(PY), "-u", str(SCRIPT)], cwd=str(ROOT), env=env,
                           capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=1800)
    except subprocess.TimeoutExpired:
        print(f"[{tag}] 超时")
        return None, False, None
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
    tf = re.findall(r"failed to converge at time : ([\d.]+)", out)
    nj = (OUT / f"{tag}.febio.log")
    njv = len(re.findall(r"[Nn]egative jacobian", nj.read_text(
        encoding="utf-8", errors="replace"))) if nj.exists() else 0
    print(f"[{tag}] NORMAL={str(normal):<5} 负J={njv:<4} 失收敛t={tf[:1] or '—'} F={F}")
    return None, normal, F


print("=== A1+A2 (SUBTALAR=1 + spring + MAX_UPS=0 + PENALTY) 量程 ===")
for d in ("+1.6", "+1.0", "+0.5", "+0.25"):
    run(f"range_A1A2_{d}", {"SUBTALAR": "1", "TALUS_SUPPORT": "spring",
                            "TALUS_K": "500", "MAX_UPS": "0", "DISP_MM": d})
print("\n=== 对照：同载荷下基线 ===")
for d in ("+1.6", "+0.5"):
    run(f"range_base_{d}", {"SUBTALAR": "0", "TALUS_SUPPORT": "fix",
                            "DISP_MM": d})
