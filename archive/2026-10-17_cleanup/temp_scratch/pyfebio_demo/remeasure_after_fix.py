"""韧带几何修复后重测 C2/C3（修正上一个脚本的 `"POSE" in kw` 键判断 bug）。"""
import os
import re
import shutil
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WORK = ROOT / "temp/thums/febio_lig"
C2D = WORK / "c2"
PY = ROOT / ".venv/Scripts/python.exe"
SCRIPT = ROOT / "scripts/ankle_fe/ankle_lig_feb.py"
STRAIN = ROOT / "scripts/ankle_fe/lig_strain.py"

BASE = {
    "DET_MIN": "0.0", "TRUNC_MM": "60", "AXIS_SRC": "full", "ALIGN_AXIS": "1",
    "FIX_TALUS": "transverse", "DRIVER": "dispx", "DISP_MM": "-1.0",
    "LOAD_SIGN": "+1", "CART_NEW": "1", "TIE": "1",
    "CART_DIR": str(ROOT / "temp/thums/cart_m1"), "LOG_DATA": "1",
}
VALS = [-3, -2, -1, 0, 1, 2, 3, 4, 5, 6]
LAT = ("ATFL", "CFL", "PTFL")
MED = ("MD1", "MD2", "MD3", "MD4")

import numpy as np


def parse_coords(p):
    t = Path(p).read_text(encoding="latin-1", errors="replace")
    return {int(m.group(1)): np.array([float(x) for x in m.group(2).split(",")])
            for m in re.finditer(r'<node id="(\d+)">([-\d.eE+,]+)</node>', t)}


def last_disp(p):
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
    if best is None:
        return None, None
    return best[0], {int(r[0]): np.array(r[1:4]) for r in best[1]}


xyz_n = parse_coords(C2D / "IE_0.feb")
print(f"{'IE':>5} {'rc':>3} {'RESULT':>7} {'|Fz|':>10} {'外侧(c)%':>9} "
      f"{'内侧(c)%':>9} {'外-内':>8}  {'CFL':>7} {'ATFL':>7} {'PTFL':>7}")
print("─" * 100)
for v in VALS:
    env = os.environ.copy()
    env.update(BASE)
    env["POSE_IE"] = str(v)
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([str(PY), "-u", str(SCRIPT)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=3600)
    out = (r.stdout or "") + (r.stderr or "")
    tag = f"lg_{v}"
    (ROOT / f"temp/pyfebio_demo/_{tag}.log").write_text(out, encoding="utf-8")
    nm = "NORMAL" if "N O R M A L   T E R M I N A T I O N" in out else "ERROR"
    fx = re.search(r"两顶面合力 = \[\s*[-\d.e+]+\s+[-\d.e+]+\s+([-\d.e+]+)", out)
    lat_s = med_s = float("nan")
    cfl = atfl = ptfl = float("nan")
    if nm == "NORMAL":
        shutil.copy2(WORK / "ankle_lig_WK.feb", C2D / f"{tag}.feb")
        shutil.copy2(WORK / "disp_WK.txt", C2D / f"disp_{tag}.txt")
        shutil.copy2(WORK / "lig_read_WK.npz", C2D / f"lig_{tag}.npz")
        xyz_p = parse_coords(C2D / f"{tag}.feb")
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
                Lc = float(np.linalg.norm((xyz_p[b] + disp[b]) - (xyz_p[a] + disp[a])))
                vals.append((Lc - L0) / L0 * 100.0)
            tot[nmx] = float(np.mean(vals)) if vals else float("nan")
        lat_s = float(np.nanmean([tot[n] for n in names if any(k in n for k in LAT)]))
        med_s = float(np.nanmean([tot[n] for n in names if any(k in n for k in MED)]))
        cfl = tot.get("CFL_calcaneofibular", float("nan"))
        atfl = tot.get("ATFL_talofibular_anterior", float("nan"))
        ptfl = tot.get("PTFL_talofibular_posterior", float("nan"))
    print(f"{v:>+5.1f} {r.returncode:>3} {nm:>7} "
          f"{(float(fx.group(1)) if fx else float('nan')):>10.2f} "
          f"{lat_s:>9.3f} {med_s:>9.3f} {lat_s-med_s:>8.3f}  "
          f"{cfl:>7.2f} {atfl:>7.2f} {ptfl:>7.2f}")
print("\n★ 对比修复前（旧韧带）: +5 外侧 +7.10 / CFL +11.79 ；外翻可用下限 −3")
