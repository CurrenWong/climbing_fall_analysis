"""A1 诊断矩阵 v2：加 MAX_UPS=0（关掉 BFGS 更新）等开关，看是否是拟牛顿更新触发的失控。"""
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

CASES = {
    "cal_noupd": {"SUBTALAR": "calcaneus", "MAX_UPS": "0"},
    "both_noupd": {"MAX_UPS": "0"},
    "cal_noupd_r2": {"SUBTALAR": "calcaneus", "MAX_UPS": "0"},
}


def run(tag, extra):
    env = os.environ.copy()
    env.update(BASE); env.update(extra); env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([str(PY), "-u", str(SCRIPT)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=1800)
    out = (r.stdout or "") + (r.stderr or "")
    (OUT / f"{tag}.stdout.txt").write_text(out, encoding="utf-8")
    m = re.search(r"FEBio rc=(-?\d+)\s+用时 ([\d.]+)s", out)
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
        rt = WORK / "reac_tibia_WK.txt"
        rf = WORK / "reac_fibula_WK.txt"
        F = sum_top(rt) + sum_top(rf)
    print(f"[{tag}] rc={rc} NORMAL={normal} 两顶面合力 z={F}")


def sum_top(p):
    import numpy as np
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


for tag, extra in CASES.items():
    run(tag, extra)
