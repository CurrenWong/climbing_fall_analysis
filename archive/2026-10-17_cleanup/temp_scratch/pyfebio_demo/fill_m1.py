"""补 IE=-1 缺格（500N）+ 落盘修正后的 C3。"""
import io
import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WORK = ROOT / "temp/thums/febio_lig"
C2D = WORK / "c2"
PY = ROOT / ".venv/Scripts/python.exe"
SCRIPT = ROOT / "scripts/ankle_fe/ankle_lig_feb.py"
W = r"C:\Users\curre\.copaw\workspaces\default"

BASE = {
    "DET_MIN": "0.0", "TRUNC_MM": "60", "AXIS_SRC": "full", "ALIGN_AXIS": "1",
    "FIX_TALUS": "transverse", "DRIVER": "dispx", "CART_NEW": "1", "TIE": "1",
    "CART_DIR": str(ROOT / "temp/thums/cart_m1"), "LOG_DATA": "1",
}


def coords(p):
    t = Path(p).read_text(encoding="latin-1", errors="replace")
    return {int(m.group(1)): np.array([float(x) for x in m.group(2).split(",")])
            for m in re.finditer(r'<node id="(\d+)">([-\d.eE+,]+)</node>', t)}


def rec(p):
    rows, tv, best = [], None, None
    for l in Path(p).read_text(encoding="utf-8", errors="replace").split("\n"):
        if l.startswith("*Step"):
            if rows and tv is not None and (best is None or abs(tv - 1) < abs(best[0] - 1)):
                best = (tv, list(rows))
            rows, tv = [], None
        elif l.startswith("*Time"):
            tv = float(l.split("=")[1])
        elif l.strip() and not l.startswith("*"):
            try:
                rows.append([float(x) for x in l.split(",")])
            except ValueError:
                pass
    if rows and tv is not None and (best is None or abs(tv - 1) < abs(best[0] - 1)):
        best = (tv, rows)
    return None if best is None else {int(r[0]): np.array(r[1:4]) for r in best[1]}


xyz_n = coords(C2D / "IE_0.feb")
for d in (1.15, 1.25, 1.35):
    env = os.environ.copy()
    env.update(BASE)
    env["POSE_IE"] = "-1"
    env["DISP_MM"] = f"+{d}"
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([str(PY), "-u", str(SCRIPT)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=3600)
    out = (r.stdout or "") + (r.stderr or "")
    (ROOT / f"temp/pyfebio_demo/_m1_{d}.log").write_text(out, encoding="utf-8")
    ok = ("FEBio rc=0" in out and "failed to converge" not in out
          and not re.search(r"negative jacobian", out, re.I))
    rr, rf = rec(WORK / "reac_tibia_WK.txt"), rec(WORK / "reac_fibula_WK.txt")
    F = (abs(sum(x[2] for x in rr.values())
             + (sum(x[2] for x in rf.values()) if rf else 0.0))) if (rr and ok) else float("nan")
    print(f"IE=-1 DISP=+{d}: rc0={ok}  F={F:.2f}")
    if ok and abs(F - 500) < 25:
        shutil.copy2(WORK / "ankle_lig_WK.feb", C2D / "_tmp.feb")
        shutil.copy2(WORK / "disp_WK.txt", C2D / "_tmp_disp.txt")
        shutil.copy2(WORK / "lig_read_WK.npz", C2D / "_tmp_lig.npz")
        xp = coords(C2D / "_tmp.feb")
        dp = rec(C2D / "_tmp_disp.txt")
        z = np.load(C2D / "_tmp_lig.npz", allow_pickle=True)
        names = [str(x) for x in z["lig_names"]]
        pr, npr = z["lig_pairs"], z["lig_npair"]
        tot = {}
        for i, nmx in enumerate(names):
            vs = []
            for j in range(int(npr[i])):
                a, b = int(pr[i, j, 0]), int(pr[i, j, 1])
                if a < 0 or b < 0:
                    continue
                L0 = float(np.linalg.norm(xyz_n[b] - xyz_n[a]))
                Lc = float(np.linalg.norm((xp[b] + dp[b]) - (xp[a] + dp[a])))
                vs.append((Lc - L0) / L0 * 100.0)
            tot[nmx] = float(np.mean(vs)) if vs else float("nan")
        LAT = ("ATFL", "CFL", "PTFL")
        MED = ("MD1", "MD2", "MD3", "MD4")
        la = float(np.nanmean([tot[n] for n in names if any(k in n for k in LAT)]))
        me = float(np.nanmean([tot[n] for n in names if any(k in n for k in MED)]))
        print(f"   ⇒ 外={la:.3f} 内={me:.3f} 外-内={la-me:.3f} "
              f"CFL={tot.get('CFL_calcaneofibular'):.2f} "
              f"ATFL={tot.get('ATFL_talofibular_anterior'):.2f}")
        break
