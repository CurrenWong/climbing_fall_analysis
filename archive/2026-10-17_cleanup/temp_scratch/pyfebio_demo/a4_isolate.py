"""A4 隔离实验：是谁把量程从 +1.6 压到 +1.2？

T1: A2-only（不接距下）@ +1.6mm，用工作配置 ⇒ 若收敛 ⇒ 距下是元凶
T2: A1+A2 @ +1.6mm，TALUS_K=5000（硬弹簧，近似恢复约束）⇒ 若收敛 ⇒ 是"约束太软"
T3: A1-only（不解除约束）@ +1.6mm + PENALTY + MAX_UPS=0 ⇒ 已知失败（对照）
T4: A2-only @ +1.6mm + TALUS_K=500（复测 E3 但用 PENALTY 口径）
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
BASE = dict(DRIVER="dispx", DISP_MM="+1.6", FIX_TALUS="transverse",
            DET_MIN="0.0", TRUNC_MM="60", AXIS_SRC="full", ALIGN_AXIS="1",
            CART_NEW="1", TIE="1", CART_DIR=str(ROOT / "temp/thums/cart_m1"),
            POSE_IE="0", LOG_DATA="1", SUBTALAR_LAUGON="PENALTY",
            MAX_UPS="0")


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
    r = subprocess.run([str(PY), "-u", str(SCRIPT)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=1800)
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
    print(f"[{tag}] NORMAL={str(normal):<5} F={F}")
    return normal, F


CASES = {
    "T1_a2only_1.6":  {"SUBTALAR": "0", "TALUS_SUPPORT": "spring",
                       "TALUS_K": "500"},
    "T2_a1a2_k5000":  {"SUBTALAR": "1", "TALUS_SUPPORT": "spring",
                       "TALUS_K": "5000"},
    "T3_a1only_1.6":  {"SUBTALAR": "1", "TALUS_SUPPORT": "fix"},
    "T4_a2only_k50":  {"SUBTALAR": "0", "TALUS_SUPPORT": "spring",
                       "TALUS_K": "50"},
}
print("=== @+1.6mm（基线零回归 = 592.6020 N）===")
for tag, extra in CASES.items():
    run(tag, extra)
