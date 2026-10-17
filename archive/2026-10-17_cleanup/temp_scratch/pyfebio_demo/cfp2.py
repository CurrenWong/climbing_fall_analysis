"""恒定载荷姿势对比 v2：修正收敛判据 + 目标 500 N（全姿势可达）。"""
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
TARGET = float(os.environ.get("TARGET_N", "500"))

BASE = {
    "DET_MIN": "0.0", "TRUNC_MM": "60", "AXIS_SRC": "full", "ALIGN_AXIS": "1",
    "FIX_TALUS": "transverse", "DRIVER": "dispx", "CART_NEW": "1", "TIE": "1",
    "CART_DIR": str(ROOT / "temp/thums/cart_m1"), "LOG_DATA": "1",
}
LAT = ("ATFL", "CFL", "PTFL")
MED = ("MD1", "MD2", "MD3", "MD4")
FAIL = {"CFL_calcaneofibular": (9.0, 13.0), "ATFL_talofibular_anterior": (11.0, 15.0)}


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


def converged(out):
    """★ 稳健判据：rc=0 且无 'failed to converge' 且无负 Jacobian。"""
    m = re.search(r"FEBio rc=(-?\d+)", out)
    return (m is not None and m.group(1) == "0"
            and "failed to converge" not in out
            and not re.search(r"negative jacobian", out, re.I))


xyz_n = coords(C2D / "IE_0.feb")
LOG = []


def run(v, d):
    env = os.environ.copy()
    env.update(BASE)
    env["POSE_IE"] = str(v)
    env["DISP_MM"] = f"{d:+.4f}"
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([str(PY), "-u", str(SCRIPT)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=3600)
    out = (r.stdout or "") + (r.stderr or "")
    LOG.append((v, d, out))
    if not converged(out):
        return None
    F = float("nan")
    rr, rf = rec(WORK / "reac_tibia_WK.txt"), rec(WORK / "reac_fibula_WK.txt")
    if rr:
        F = abs(sum(x[2] for x in rr.values())
                + (sum(x[2] for x in rf.values()) if rf else 0.0))
    shutil.copy2(WORK / "ankle_lig_WK.feb", C2D / "_tmp.feb")
    shutil.copy2(WORK / "disp_WK.txt", C2D / "_tmp_disp.txt")
    shutil.copy2(WORK / "lig_read_WK.npz", C2D / "_tmp_lig.npz")
    return F


def strain_at():
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
    return tot, names


print(f"=== 恒定载荷 {TARGET:.0f} N 下的姿势对比（割线法）===")
print(f"{'IE':>5} {'DISP':>7} {'力(N)':>8} {'外(c)%':>8} {'内(c)%':>8} {'外-内':>8} "
      f"{'CFL':>7} {'CFL风险':>9} {'ATFL':>7} {'ATFL风险':>10} {'迭代':>4}")
print("─" * 104)
for v in (-3, -2, -1, 0, 1, 2, 3, 4, 5):
    d, hist, done = 1.5, [], False
    F = None
    for it in range(7):
        F = run(v, d)
        if F is None or F != F:
            # 失败 ⇒ 回退步长
            d = max(0.2, d * 0.7)
            continue
        hist.append((d, F))
        if abs(F - TARGET) / TARGET < 0.015:
            done = True
            break
        if len(hist) >= 2 and hist[-1][1] != hist[-2][1]:
            slope = (hist[-1][1] - hist[-2][1]) / (hist[-1][0] - hist[-2][0])
            d = hist[-1][0] + (TARGET - hist[-1][1]) / slope
        else:
            d = hist[-1][0] * TARGET / max(hist[-1][1], 1e-9)
        d = min(max(d, 0.05), 4.0)
    if F is None or F != F:
        print(f"{v:>+5.0f}   —— 未能收敛到目标（力提取失败）")
        continue
    tot, names = strain_at()
    la = float(np.nanmean([tot[n] for n in names if any(k in n for k in LAT)]))
    me = float(np.nanmean([tot[n] for n in names if any(k in n for k in MED)]))
    c, a = tot.get("CFL_calcaneofibular", np.nan), tot.get("ATFL_talofibular_anterior", np.nan)
    rl, rh = FAIL["CFL_calcaneofibular"]
    al, ah = FAIL["ATFL_talofibular_anterior"]
    print(f"{v:>+5.0f} {d:>+7.3f} {F:>8.2f} {la:>8.3f} {me:>8.3f} {la-me:>8.3f} "
          f"{c:>7.2f} {c/rh:>4.2f}/{c/rl:<4.2f} {a:>7.2f} {a/ah:>4.2f}/{a/al:<4.2f} "
          f"{len(hist):>4}{'' if done else ' ⚠'}")
print(f"\n★ 风险比 = 总应变/失效应变（中值/保守下限）。≥1 ⇒ 越阈。")
