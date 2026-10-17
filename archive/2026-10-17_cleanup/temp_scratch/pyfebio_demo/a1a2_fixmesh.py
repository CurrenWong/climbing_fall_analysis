"""关键假设测试：用 cart_m1_fix（已修掉 11 个翻转单元的网格）跑 A1+A2。

E2p_fix : SUBTALAR=1 + PENALTY + spring + MAX_UPS=0 + CART_DIR=cart_m1_fix
E1p_fix : SUBTALAR=1 + PENALTY + fix + MAX_UPS=0 + CART_DIR=cart_m1_fix
E0_fix  : 基线配置 + CART_DIR=cart_m1_fix（对照：修网格对基线的影响）
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
FIX = str(ROOT / "temp/thums/cart_m1_fix")
OLD = str(ROOT / "temp/thums/cart_m1")
BASE = dict(DRIVER="dispx", DISP_MM="+1.6", DET_MIN="0.0", TRUNC_MM="60",
            AXIS_SRC="full", ALIGN_AXIS="1", FIX_TALUS="transverse",
            CART_NEW="1", TIE="1", POSE_IE="0", LOG_DATA="1",
            SUBTALAR_LAUGON="PENALTY")


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
    # 负 J 的两种写法都要抓
    nj = len(re.findall(r"[Nn]egative jacobian", (OUT / f"{tag}.febio.log")
                        .read_text(encoding="utf-8", errors="replace")
                        if (OUT / f"{tag}.febio.log").exists() else ""))
    tf = re.findall(r"failed to converge at time : ([\d.]+)", out)
    print(f"[{tag}] rc={rc} NORMAL={normal} 负J={nj} 失收敛t={tf[:1] or '—'} F={F}")
    return rc, normal, F


CASES = {
    "E2p_fix": {"SUBTALAR": "1", "TALUS_SUPPORT": "spring", "TALUS_K": "500",
                "MAX_UPS": "0", "CART_DIR": FIX},
    "E1p_fix": {"SUBTALAR": "1", "TALUS_SUPPORT": "fix",
                "MAX_UPS": "0", "CART_DIR": FIX},
    "E0_fix":  {"SUBTALAR": "0", "TALUS_SUPPORT": "fix", "CART_DIR": FIX},
}
res = {}
for tag, extra in CASES.items():
    res[tag] = run(tag, extra)
print("\n=== 汇总（对照：E0=592.6020 / E1p 0.645 / E2p 0.761）===")
for k, v in res.items():
    print(f"  {k:<10} rc={v[0]} NORMAL={v[1]} F={v[2]}")
