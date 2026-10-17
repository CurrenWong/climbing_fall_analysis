"""A1+A2 续推：停滞 ⇒ 试 TIME_STEPS↑ / MAX_REFS↑（E1 与 E3 各试两档）。"""
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
BASE = dict(DRIVER="dispx", DISP_MM="+1.6", DET_MIN="0.0", TRUNC_MM="60",
            AXIS_SRC="full", ALIGN_AXIS="1", FIX_TALUS="transverse",
            CART_NEW="1", TIE="1", CART_DIR="temp/thums/cart_m1",
            POSE_IE="0", LOG_DATA="1", MAX_UPS="0")


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
    t_fail = re.findall(r"failed to converge at time : ([\d.]+)", out)
    print(f"[{tag}] rc={rc} NORMAL={normal} 失收敛t={t_fal if (t_fal:=t_fail[:1]) else '—'} F={F}")
    return rc, normal, F


CASES = {
    "E1_ts40":  {"SUBTALAR": "1", "TALUS_SUPPORT": "fix",
                 "TIME_STEPS": "40", "MAX_REFS": "100"},
    "E1_ts40r150": {"SUBTALAR": "1", "TALUS_SUPPORT": "fix",
                    "TIME_STEPS": "40", "MAX_REFS": "150"},
    "E3_ts40":  {"SUBTALAR": "0", "TALUS_SUPPORT": "spring",
                 "TALUS_K": "500", "TIME_STEPS": "40", "MAX_REFS": "100"},
}
res = {}
for tag, extra in CASES.items():
    res[tag] = run(tag, extra)
print("\n=== 汇总 ===")
for k, v in res.items():
    print(f"  {k:<14} rc={v[0]} NORMAL={v[1]} F={v[2]}")
