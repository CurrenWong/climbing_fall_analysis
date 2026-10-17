"""恒定载荷口径：对每个姿势二分/割线求 DISP 使关节力 = 600 N，再比较应变。

意义：固定位移下比较会**混淆载荷与姿势**（姿势本身让关节变软 ⇒ 力下降）。
用户目标（判不同姿势的风险）要求**同一载荷**下比较。
"""
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
TARGET = float(os.environ.get("TARGET_N", "600"))

BASE = {
    "DET_MIN": "0.0", "TRUNC_MM": "60", "AXIS_SRC": "full", "ALIGN_AXIS": "1",
    "FIX_TALUS": "transverse", "DRIVER": "dispx", "CART_NEW": "1", "TIE": "1",
    "CART_DIR": str(ROOT / "temp/thums/cart_m1"), "LOG_DATA": "1",
}
LAT = ("ATFL", "CFL", "PTFL")
MED = ("MD1", "MD2", "MD3", "MD4")
FAIL = {"CFL_calcaneofibular": (9.0, 13.0), "ATFL_talofibular_anterior": (11.0, 15.0),
        "PTFL_talofibular_posterior": (None, None)}


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
    ok = "N O R M A L   T E R M I N A T I O N" in out
    if not ok:
        return None, None, out
    F = float("nan")
    rr = rec(WORK / "reac_tibia_WK.txt")
    rf = rec(WORK / "reac_fibula_WK.txt")
    if rr:
        F = abs(sum(x[2] for x in rr.values()) + (sum(x[2] for x in rf.values())
                                                 if rf else 0.0))
    shutil.copy2(WORK / "ankle_lig_WK.feb", C2D / "_tmp.feb")
    shutil.copy2(WORK / "disp_WK.txt", C2D / "_tmp_disp.txt")
    shutil.copy2(WORK / "lig_read_WK.npz", C2D / "_tmp_lig.npz")
    return F, out, out


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


print(f"=== 恒定载荷 {TARGET:.0f} N 下的姿势对比（割线法求 DISP）===")
print(f"{'IE':>5} {'DISP':>7} {'力(N)':>9} {'外(c)%':>8} {'内(c)%':>8} {'外-内':>8} "
      f"{'CFL':>7} {'CFL风险比':>10} {'ATFL':>7} {'ATFL风险比':>11}")
print("─" * 104)
for v in (-3, 0, 3, 5):
    d, F = 2.0, None
    hist = []
    for it in range(6):
        F, _, out = run(v, d)
        if F is None or F != F:
            print(f"{v:>+5.0f}  迭代 {it}: DISP={d:+.3f} 失败")
            break
        hist.append((d, F))
        if abs(F - TARGET) / TARGET < 0.01:
            break
        # 割线
        if len(hist) >= 2 and hist[-1][1] != hist[-2][1]:
            slope = (hist[-1][1] - hist[-2][1]) / (hist[-1][0] - hist[-2][0])
            d = hist[-1][0] + (TARGET - hist[-1][1]) / slope
        else:
            d = hist[-1][0] * TARGET / hist[-1][1]
        d = min(max(d, 0.1), 6.0)
    if F is None or F != F:
        continue
    tot, names = strain_at()
    la = float(np.nanmean([tot[n] for n in names if any(k in n for k in LAT)]))
    me = float(np.nanmean([tot[n] for n in names if any(k in n for k in MED)]))
    c = tot.get("CFL_calcaneofibular", float("nan"))
    a = tot.get("ATFL_talofibular_anterior", float("nan"))
    rl, rh = FAIL["CFL_calcaneofibular"]
    cr = f"{c/rh:.2f}/{c/rl:.2f}" if rl else "—"
    al, ah = FAIL["ATFL_talofibular_anterior"]
    ar = f"{a/ah:.2f}/{a/al:.2f}" if al else "—"
    print(f"{v:>+5.0f} {d:>+7.3f} {F:>9.2f} {la:>8.3f} {me:>8.3f} {la-me:>8.3f} "
          f"{c:>7.2f} {cr:>10} {a:>7.2f} {ar:>11}")
print("\n★ 风险比 = 总应变/失效应变（中值/保守下限）。≥1 ⇒ 越阈。")
print("★ 这才是「同一载荷下不同姿势」的干净比较。")
