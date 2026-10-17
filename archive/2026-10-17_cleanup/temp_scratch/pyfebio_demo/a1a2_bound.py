"""A1+A2 收敛边界扫描：找能到达的最大 DISP（600 N 目标对应 ~+1.3~1.4mm）。"""
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
BASE = dict(DRIVER="dispx", FIX_TALUS="transverse", DET_MIN="0.0", TRUNC_MM="60",
            AXIS_SRC="full", ALIGN_AXIS="1", CART_NEW="1", TIE="1",
            CART_DIR=str(ROOT / "temp/thums/cart_m1"), POSE_IE="0",
            LOG_DATA="1", SUBTALAR_LAUGEN="PENALTY", SUBTALAR_LAUGON="PENALTY",
            SUBTALAR="1", TALUS_SUPPORT="spring", TALUS_K="500", MAX_UPS="0")


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


rows = []
for d in ("+1.1", "+1.2", "+1.3", "+1.35", "+1.4", "+1.5"):
    nm, F = run(f"bnd_{d}", {"DISP_MM": d})
    rows.append((d, nm, F))
print("\n" + "=" * 70)
print("A1+A2 收敛边界（DISP_MM →）")
print("=" * 70)
ok = [r for r in rows if r[1]]
for d, nm, F in rows:
    print(f"  {d:>6}  NORMAL={str(nm):<5} F={F}")
if ok:
    dmax = max(ok, key=lambda r: float(r[0]))[0]
    Fmax = max(ok, key=lambda r: float(r[0]))[2]
    print(f"\n  ⇒ 收敛上界 ≈ {dmax} mm，对应 {Fmax:.2f} N"
          f" ⇒ {'✅ 可覆盖 600 N' if Fmax >= 600 else '❌ 够不到 600 N'}")
