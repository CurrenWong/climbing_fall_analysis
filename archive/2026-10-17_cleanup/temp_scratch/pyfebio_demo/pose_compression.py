"""正确符号（压缩 DISP_MM>0）下重测姿势矩阵的韧带应变剂量-响应。"""
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

DISP = os.environ.get("CDISP", "+1.0")          # ★ 压缩方向
BASE = {
    "DET_MIN": "0.0", "TRUNC_MM": "60", "AXIS_SRC": "full", "ALIGN_AXIS": "1",
    "FIX_TALUS": "transverse", "DRIVER": "dispx", "DISP_MM": DISP,
    "CART_NEW": "1", "TIE": "1",
    "CART_DIR": str(ROOT / "temp/thums/cart_m1"), "LOG_DATA": "1",
}
VALS = [-3, -2, -1, 0, 1, 2, 3, 4, 5]
LAT = ("ATFL", "CFL", "PTFL")
MED = ("MD1", "MD2", "MD3", "MD4")


def coords(p):
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
        return None
    return {int(r[0]): np.array(r[1:4]) for r in best[1]}


NEU = C2D / "IE_0.feb"
if not NEU.exists():
    print(f"✗ 中性参照缺失 {NEU}（先跑 check_c2.py）")
    raise SystemExit(1)
xyz_n = coords(NEU)
print(f"中性参照 {NEU.name} {len(xyz_n)} 节点   压缩 DISP_MM={DISP}\n")
print(f"{'IE':>5} {'rc':>3} {'RESULT':>7} {'|反力z|':>9} {'外(c)%':>8} "
      f"{'内(c)%':>8} {'外-内':>8} {'CFL':>7} {'ATFL':>7} {'PTFL':>7}")
print("─" * 96)
for v in VALS:
    env = os.environ.copy()
    env.update(BASE)
    env["POSE_IE"] = str(v)
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([str(PY), "-u", str(SCRIPT)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=3600)
    out = (r.stdout or "") + (r.stderr or "")
    tg = f"cp_{v}"
    (ROOT / f"temp/pyfebio_demo/_{tg}.log").write_text(out, encoding="utf-8")
    nm = "NORMAL" if "N O R M A L   T E R M I N A T I O N" in out else "ERROR"
    nj = len(re.findall(r"negative jacobian", out, re.I))
    la = me = cfl = atfl = ptfl = float("nan")
    Fz = float("nan")
    if nm == "NORMAL":
        shutil.copy2(WORK / "ankle_lig_WK.feb", C2D / f"{tg}.feb")
        shutil.copy2(WORK / "disp_WK.txt", C2D / f"disp_{tg}.txt")
        shutil.copy2(WORK / "lig_read_WK.npz", C2D / f"lig_{tg}.npz")
        xp = coords(C2D / f"{tg}.feb")
        dp = last_disp(C2D / f"disp_{tg}.txt")
        z = np.load(C2D / f"lig_{tg}.npz", allow_pickle=True)
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
        la = float(np.nanmean([tot[n] for n in names if any(k in n for k in LAT)]))
        me = float(np.nanmean([tot[n] for n in names if any(k in n for k in MED)]))
        cfl = tot.get("CFL_calcaneofibular", float("nan"))
        atfl = tot.get("ATFL_talofibular_anterior", float("nan"))
        ptfl = tot.get("PTFL_talofibular_posterior", float("nan"))
        rp = WORK / "reac_tibia_WK.txt"
        if rp.exists():
            rr = last_disp(rp)
            if rr:
                Fz = float(abs(sum(x[2] for x in rr.values())))
    print(f"{v:>+5.1f} {r.returncode:>3} {nm:>7} {Fz:>9.2f} {la:>8.3f} "
          f"{me:>8.3f} {la-me:>8.3f} {cfl:>7.2f} {atfl:>7.2f} {ptfl:>7.2f}"
          + (f"   ⚠{nj}负J" if nj else ""))
print("\n★ 对比（旧·**拉伸**方向 DISP_MM=-1.0）：+5 外-内=+9.07/CFL=+11.79；"
      "-3 外-内=-6.76")
