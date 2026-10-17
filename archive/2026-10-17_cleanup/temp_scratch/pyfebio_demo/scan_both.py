"""双侧姿势可用范围扫描 + 两侧韧带应变（(c) 总应变，供 C3 用）。"""
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

BASE = {
    "DET_MIN": "0.0", "TRUNC_MM": "60", "AXIS_SRC": "full", "ALIGN_AXIS": "1",
    "FIX_TALUS": "transverse", "DRIVER": "dispx", "DISP_MM": "-1.0",
    "LOAD_SIGN": "+1", "CART_NEW": "1", "TIE": "1",
    "CART_DIR": str(ROOT / "temp/thums/cart_m1"), "LOG_DATA": "1",
}
VALS = [-5, -4, -3.5, -3, -2, -1, 1, 2, 3, 4, 5, 6, 8]
LAT = ("ATFL", "CFL", "PTFL")
MED = ("MD1", "MD2", "MD3", "MD4")


def parse_feb_coords(p):
    t = p.read_text(encoding="latin-1", errors="replace")
    out = {}
    for m in re.finditer(r'<node id="(\d+)">([-\d.eE+,]+)</node>', t):
        out[int(m.group(1))] = np.array(
            [float(x) for x in m.group(2).split(",")], float)
    return out


def last_disp(p):
    rows, tv, best = [], None, None
    for l in p.read_text(encoding="utf-8", errors="replace").split("\n"):
        if l.startswith("*Step"):
            if rows and tv is not None:
                if best is None or abs(tv - 1) < abs(best[0] - 1):
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
    if best is None:
        return None, None
    return best[0], {int(r[0]): np.array(r[1:4]) for r in best[1]}


xyz_n = parse_feb_coords(C2D / "IE_0.feb")

print(f"{'IE':>5} {'rc':>3} {'NORMAL':>7} {'负J':>4} {'|Fz|':>9} "
      f"{'外侧(c)%':>9} {'内侧(c)%':>9} {'外-内':>8}  {'外侧韧带明细':<34}")
print("─" * 118)
rows = []
for v in VALS:
    tag = f"IE_{int(round(v*100))}"
    env = os.environ.copy()
    env.update(BASE)
    env["POSE_IE"] = str(v)
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([str(PY), "-u", str(SCRIPT)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=3600)
    out = (r.stdout or "") + (r.stderr or "")
    io.open(ROOT / f"temp/pyfebio_demo/_scan_{tag}.log", "w",
            encoding="utf-8").write(out)
    nm = "NORMAL" if "N O R M A L   T E R M I N A T I O N" in out else "ERROR"
    nj = len(re.findall(r"negative jacobian", out, re.I))
    fz = re.search(r"两顶面合力 = \[\s*[-\d.e+]+\s+[-\d.e+]+\s+([-\d.e+]+)", out)
    lat_s = med_s = float("nan")
    detail = ""
    if nm == "NORMAL":
        pf = C2D / f"{tag}.feb"
        shutil.copy2(WORK / "ankle_lig_WK.feb", pf)
        shutil.copy2(WORK / "disp_WK.txt", C2D / f"disp_{tag}.txt")
        shutil.copy2(WORK / "lig_read_WK.npz", C2D / f"lig_{tag}.npz")
        xyz_p = parse_feb_coords(pf)
        t0, disp = last_disp(C2D / f"disp_{tag}.txt")
        z = np.load(C2D / f"lig_{tag}.npz", allow_pickle=True)
        names = [str(x) for x in z["lig_names"]]
        pairs, npair = z["lig_pairs"], z["lig_npair"]
        tot = {}
        for i, nmx in enumerate(names):
            vals = []
            for j in range(int(npair[i])):
                a, b = int(pairs[i, j, 0]), int(pairs[i, j, 1])
                if a < 0 or b < 0:
                    continue
                L0 = float(np.linalg.norm(xyz_n[b] - xyz_n[a]))
                Lc = float(np.linalg.norm(
                    (xyz_p[b] + disp[b]) - (xyz_p[a] + disp[a])))
                vals.append((Lc - L0) / L0 * 100.0)
            tot[nmx] = float(np.mean(vals)) if vals else float("nan")
        lat_s = float(np.nanmean([tot[n] for n in names if any(k in n for k in LAT)]))
        med_s = float(np.nanmean([tot[n] for n in names if any(k in n for k in MED)]))
        detail = " ".join(f"{n.split('_')[0]}:{tot[n]:+.2f}"
                          for n in names if any(k in n for k in LAT))
    print(f"{v:>+5.1f} {r.returncode:>3} {nm:>7} {nj:>4} "
          f"{(float(fz.group(1)) if fz else float('nan')):>9.2f} "
          f"{lat_s:>9.3f} {med_s:>9.3f} {lat_s-med_s:>8.3f}  {detail:<34}")
    rows.append((v, nm, lat_s, med_s))

print("\n★ 判据（③）：内翻侧 (外−内) 应显著为正；外翻侧应显著为负。")
print("★ 阈值（C3 待定）：用 (c) 总应变，结合解剖失效应变（ATFL 11~15%、CFL 9~13%）。")
