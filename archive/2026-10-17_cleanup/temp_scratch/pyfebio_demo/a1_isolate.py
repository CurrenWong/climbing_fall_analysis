"""A1 起效量的隔离实验（正对照）：SUBTALAR=0 + MAX_UPS=0。

目的：把「接入距下关节」的效应与「MAX_UPS=0 求解器设置」的效应分开。
  F(不接, BFGS 开) = 592.6020   （已知基线）
  F(不接, BFGS 关) = ?          ← 本脚本测
  F(接入, BFGS 关) = 592.5697   （a1_matrix2.py 已测）
⇒ 接入的真实效应 = F(接入,关) − F(不接,关)
"""
import os
import re
import shutil
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WORK = ROOT / "temp/thums/febio_lig"
OUT = ROOT / "temp/pyfebio_demo/a1_matrix"
PY = ROOT / ".venv/Scripts/python.exe"
SCRIPT = ROOT / "scripts/ankle_fe/ankle_lig_feb.py"
BASE = dict(DRIVER="dispx", DISP_MM="+1.6", DET_MIN="0.0", TRUNC_MM="60",
            AXIS_SRC="full", ALIGN_AXIS="1", FIX_TALUS="transverse",
            CART_NEW="1", TIE="1", CART_DIR="temp/thums/cart_m1",
            POSE_IE="0", LOG_DATA="1")
import numpy as np


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
    m = re.search(r"FEBio rc=(-?\d+)", out)
    rc = int(m.group(1)) if m else None
    lg = WORK / "log_WK.txt"
    normal = lg.exists() and "N O R M A L   T E R M I N A T I O N" in \
        lg.read_text(encoding="utf-8", errors="replace")
    for src, dst in ((WORK / "febio_WK.log", f"{tag}.febio.log"),
                     (WORK / "ankle_lig_WK.feb", f"{tag}.feb")):
        if src.exists():
            shutil.copy(src, OUT / dst)
    F = None
    if normal:
        F = (sum_top(WORK / "reac_tibia_WK.txt")
             + sum_top(WORK / "reac_fibula_WK.txt"))
    print(f"[{tag}] rc={rc} NORMAL={normal} 两顶面合力 z={F}")
    return F


print("=== 隔离矩阵（每格 = 1 次求解）===")
cases = [("nosub_noupd", {"SUBTALAR": "0", "MAX_UPS": "0"}),
         ("sub_noupd_ctl", {"MAX_UPS": "0"})]
res = {}
for tag, extra in cases:
    res[tag] = run(tag, extra)

F0, F1 = 592.6020, res.get("nosub_noupd")
F2 = res.get("sub_noupd_ctl")
print("\n" + "=" * 72)
print(f"  F(不接, BFGS 开) = {F0:9.4f}   ← 已知基线")
if F1:
    print(f"  F(不接, BFGS 关) = {F1:9.4f}   ← MAX_UPS=0 自身效应 "
          f"({(F1-F0)/F0*100:+.4f}%)")
if F2:
    print(f"  F(接入, BFGS 关) = {F2:9.4f}")
if F1 and F2:
    print(f"\n⭐ **接入距下关节的真实效应** = {F2-F1:+.4f} N "
          f"= {(F2-F1)/F1*100:+.4f}%")
    print(f"   对比「缺距下关节 44.8%」的说法 ⇒ "
          f"{'❌ 被直接测量否证（差 4 个量级）' if abs((F2-F1)/F1) < 0.01 else '⚠ 需再查'}")
