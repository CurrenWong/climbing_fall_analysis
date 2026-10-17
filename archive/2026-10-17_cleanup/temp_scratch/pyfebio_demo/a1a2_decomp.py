"""A1/A2 正确分解（全部在 DISP=+1.0，四格均应收敛）。

A0 基线      : SUBTALAR=0 TALUS_SUPPORT=fix    （BFGS 开，零回归口径）
A1 仅接关节  : SUBTALAR=1 TALUS_SUPPORT=fix    MAX_UPS=0
A2 仅解除约束: SUBTALAR=0 TALUS_SUPPORT=spring MAX_UPS=0
★ A1A2       : SUBTALAR=1 TALUS_SUPPORT=spring MAX_UPS=0
A0n 基线+关BFGS: SUBTALAR=0 TALUS_SUPPORT=fix   MAX_UPS=0  ← 隔离 MAX_UPS 自身
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
BASE = dict(DRIVER="dispx", DISP_MM="+1.0", FIX_TALUS="transverse",
            DET_MIN="0.0", TRUNC_MM="60", AXIS_SRC="full", ALIGN_AXIS="1",
            CART_NEW="1", TIE="1", CART_DIR=str(ROOT / "temp/thums/cart_m1"),
            POSE_IE="0", LOG_DATA="1", SUBTALAR_LAUGON="PENALTY")


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
    "D0_base_bfgs":   {"SUBTALAR": "0", "TALUS_SUPPORT": "fix"},
    "D0n_base_noBFGS": {"SUBTALAR": "0", "TALUS_SUPPORT": "fix", "MAX_UPS": "0"},
    "D1_a1_only":     {"SUBTALAR": "1", "TALUS_SUPPORT": "fix", "MAX_UPS": "0"},
    "D2_a2_only":     {"SUBTALAR": "0", "TALUS_SUPPORT": "spring",
                       "TALUS_K": "500", "MAX_UPS": "0"},
    "D12_A1A2":       {"SUBTALAR": "1", "TALUS_SUPPORT": "spring",
                       "TALUS_K": "500", "MAX_UPS": "0"},
}
res = {}
for tag, extra in CASES.items():
    res[tag] = run(tag, extra)

print("\n" + "=" * 76)
print("DISP=+1.0 mm 分解（基线零回归值 = 418.10 N）")
print("=" * 76)
g = {k: v[1] for k, v in res.items()}
for k, (nm, F) in res.items():
    print(f"  {k:<18} NORMAL={str(nm):<5} F={F}")
if g.get("D0n_base_noBFGS") and g.get("D12_A1A2") and g.get("D0_base_bfgs"):
    b, n, x = g["D0_base_bfgs"], g["D0n_base_noBFGS"], g["D12_A1A2"]
    print(f"\n  MAX_UPS=0 自身效应      : {n-b:+.3f} N ({(n-b)/b*100:+.3f}%)")
    if g.get("D1_a1_only"):
        print(f"  仅接关节 A1 效应        : {g['D1_a1_only']-n:+.3f} N "
              f"({(g['D1_a1_only']-n)/n*100:+.3f}%)")
    if g.get("D2_a2_only"):
        print(f"  仅解除约束 A2 效应      : {g['D2_a2_only']-n:+.3f} N "
              f"({(g['D2_a2_only']-n)/n*100:+.3f}%)")
    print(f"  ★ A1+A2 合计            : {x-n:+.3f} N ({(x-n)/n*100:+.3f}%)")
