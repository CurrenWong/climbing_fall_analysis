"""A1+A2 合并实验矩阵（4 格，每格 1 次求解）。

E0 基线        : SUBTALAR=0  TALUS_SUPPORT=fix     MAX_UPS=10  → 592.6020（已知）
E1 仅 A1       : SUBTALAR=1  TALUS_SUPPORT=fix     MAX_UPS=0   → 592.5697（已知）
★ E2 A1+A2 合并 : SUBTALAR=1  TALUS_SUPPORT=spring  MAX_UPS=0
E3 仅 A2(对照)  : SUBTALAR=0  TALUS_SUPPORT=spring  MAX_UPS=0

测量：rc / NORMAL / 两顶面合力 / （若 LOG_CONTACT 可读）距下接触力。
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
OUT.mkdir(parents=True, exist_ok=True)
PY = ROOT / ".venv/Scripts/python.exe"
SCRIPT = ROOT / "scripts/ankle_fe/ankle_lig_feb.py"
BASE = dict(DRIVER="dispx", DISP_MM="+1.6", DET_MIN="0.0", TRUNC_MM="60",
            AXIS_SRC="full", ALIGN_AXIS="1", FIX_TALUS="transverse",
            CART_NEW="1", TIE="1", CART_DIR="temp/thums/cart_m1",
            POSE_IE="0", LOG_DATA="1")


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
    # 关键：从 stdout 抓距下接触/弹簧信息
    keys = []
    for ln in out.split("\n"):
        if ("subtalar" in ln.lower() and ("接" in ln or "pair" in ln.lower()
                                          or "contact" in ln.lower())) \
                or "talus-spring" in ln:
            keys.append(ln.strip()[:150])
    print(f"[{tag}] rc={rc} NORMAL={normal} F={F}")
    for k in keys[:6]:
        print(f"        {k}")
    return rc, normal, F


CASES = {
    "E0_base":        {"SUBTALAR": "0", "TALUS_SUPPORT": "fix"},
    "E1_a1_only":     {"SUBTALAR": "1", "TALUS_SUPPORT": "fix", "MAX_UPS": "0"},
    "E2_A1A2":        {"SUBTALAR": "1", "TALUS_SUPPORT": "spring",
                       "TALUS_K": "500", "MAX_UPS": "0"},
    "E3_a2_only":     {"SUBTALAR": "0", "TALUS_SUPPORT": "spring",
                       "TALUS_K": "500", "MAX_UPS": "0"},
}
res = {}
for tag, extra in CASES.items():
    res[tag] = run(tag, extra)

print("\n" + "=" * 74)
print("A1+A2 合并矩阵")
print("=" * 74)
for tag, (rc, nm, F) in res.items():
    print(f"  {tag:<14} rc={rc}  NORMAL={str(nm):<5}  F={F}")
if res.get("E2_A1A2") and res["E2_A1A2"][2] and res.get("E0_base") and res["E0_base"][2]:
    f0 = res["E0_base"][2]
    f2 = res["E2_A1A2"][2]
    print(f"\n  ★ A1+A2 合并相对基线：{f2-f0:+.4f} N = {(f2-f0)/f0*100:+.2f}%")
