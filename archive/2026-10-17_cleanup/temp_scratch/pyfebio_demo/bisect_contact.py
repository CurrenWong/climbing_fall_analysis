"""二分：关掉接触，看 tet4 与 penta6 是否还差 2.2 倍。

- 若关接触后两者一致 ⇒ 差异在**接触面**（面数/形状不同）
- 若仍差 2.2 倍 ⇒ 差异在实体网格本身
"""
import os
import re
import subprocess
from pathlib import Path
import numpy as np
import h5py
import pyfebio.xplt as X

ROOT = Path(r"D:\Project\climbing_fall_analysis")
W = ROOT / "temp/thums/febio_lig"
PY = ROOT / ".venv/Scripts/python.exe"
MODEL = ROOT / "scripts/ankle_fe/ankle_lig_feb.py"
AX = np.array([0.8153, 0.135, 0.5631])


def run(elem, no_contact):
    env = dict(os.environ)
    env.update(dict(MATSET="WK", LIG_MODE="off", DRIVER="disp",
                    TOPFIX="lateral", TRUNC_MM="60", FIX_TALUS="all",
                    TIME_STEPS="20", DISP_MM="0.1472", CART_ELEM=elem,
                    NO_CONTACT="1" if no_contact else "0"))
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([str(PY), "-u", str(MODEL)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=3600)
    o = (r.stdout or "") + (r.stderr or "")
    m = re.search(r"FEBio rc=(\d+)", o)
    rc = int(m.group(1)) if m else None
    if rc != 0:
        return rc, float("nan"), float("nan"), o
    xp = sorted(W.glob("ankle_lig_WK*.xplt"),
                key=lambda p: -p.stat().st_mtime)[0]
    X.to_hdf5(str(xp), str(W / "_nc.h5"))
    with h5py.File(W / "_nc.h5", "r") as h:
        st = sorted(int(k) for k in h["states"])
        last = max(st)
        R = np.asarray(
            h[f"states/{last}/node_data/reaction forces/all"][()],
            float).reshape(-1, 3)
        D = np.asarray(
            h[f"states/{last}/node_data/displacement/all"][()],
            float).reshape(-1, 3)
    return rc, float((R @ AX).sum()), float(np.linalg.norm(D, axis=1).max()), o


print("全部关掉韧带（LIG_MODE=off），只比实体网格；接触开/关各一轮")
print(f"{'单元':>8} {'接触':>6} {'rc':>3} {'F(N)':>10} {'|u|max':>9}")
print("-" * 46)
res = {}
for nc in (True, False):
    for elem in ("tet4", "prism"):
        rc, F, um, o = run(elem, nc)
        res[(elem, nc)] = F
        note = "关" if nc else "开"
        print(f"{elem:>8} {note:>6} {rc:>3} {F:10.1f} {um:9.4f}")
    a, b = res[("tet4", nc)], res[("prism", nc)]
    if b and b == b and b != 0:
        print(f"         -> tet4/penta6 = {a / b:.3f}   "
              f"({'两者一致' if abs(a / b - 1) < 0.15 else '仍有差异'})")
    print()
