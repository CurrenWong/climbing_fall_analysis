"""C5：韧带应变的排序对罚因子（PENALTY 0.1 vs 1）是否鲁棒。"""
import io
import os
import re
import subprocess
from pathlib import Path

import numpy as np

R = Path(r"D:\Project\climbing_fall_analysis")
W = R / "temp" / "thums" / "febio_lig"
FEB = r"D:\Program\FEBioStudio\bin\febio4.exe"
PY = str(R / ".venv" / "Scripts" / "python.exe")
COMMON = dict(DRIVER="dispx", DISP_MM="-1.0", LOAD_SIGN="+1",
              FIX_TALUS="transverse", TRUNC_MM="60", AXIS_SRC="full",
              ALIGN_AXIS="1", CART_NEW="1", TIE="1", DET_MIN="0.0",
              TALUS_RIGID="0", TAL_ROTLOCK="0", TIME_STEPS="20",
              LOG_DATA="1", LAUGON="PENALTY",
              CART_DIR=str(R / "temp" / "thums" / "cart_m1"))
LAT = ("CFL_calcaneofibular", "ATFL_talofibular_anterior",
       "PTFL_talofibular_posterior")
MED = ("MD2_tibiocalcaneal", "MD3_tibionavicular", "MD4_tibiotalar_ant",
       "MD1_tibiotalar_post")


def coords_of(fp):
    t = fp.read_text(encoding="latin-1", errors="replace")
    return {int(m.group(1)): [float(m.group(2)), float(m.group(3)),
                              float(m.group(4))]
            for m in re.finditer(r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,'
                                 r'\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)', t)}


def disp_last(fp):
    blocks = []
    for ln in io.open(fp, encoding="utf-8", errors="replace"):
        ln = ln.strip()
        if ln.startswith("*Step"):
            blocks.append([])
            continue
        if not ln or ln.startswith("*"):
            continue
        try:
            v = [float(x) for x in ln.replace(",", " ").split()]
        except ValueError:
            continue
        if blocks:
            blocks[-1].append(v)
    return {int(v[0]): np.array(v[1:4]) for v in (blocks[-1] if blocks else [])}


def measure(ie, pen):
    feb = W / "ankle_lig_WK.feb"
    if feb.exists():
        feb.unlink()
    e = dict(os.environ)
    e.update(COMMON)
    e["POSE_IE"] = ie
    e["PENALTY"] = pen
    subprocess.run([PY, "-u", str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
                   capture_output=True, text=True, env=e, cwd=str(R), timeout=900)
    if not feb.exists():
        return None
    subprocess.run([FEB, "-i", feb.name, "-o", "s.xplt"], capture_output=True,
                   text=True, cwd=str(W), timeout=1800)
    X, U = coords_of(feb), disp_last(W / "disp_WK.txt")
    z = np.load(W / "lig_read_WK.npz", allow_pickle=True)
    names = [str(x) for x in z["lig_names"]]
    PR, npr = z["lig_pairs"], z["lig_npair"]
    out = {}
    for i, nm in enumerate(names):
        dl = []
        for a, b in PR[i, :int(npr[i]), :]:
            a, b = int(a), int(b)
            if a < 0 or b < 0 or a not in X or b not in X:
                continue
            p0 = np.array(X[a]) - np.array(X[b])
            L0 = np.linalg.norm(p0)
            if L0 < 1e-9:
                continue
            p1 = p0 + (U.get(a, np.zeros(3)) - U.get(b, np.zeros(3)))
            dl.append((np.linalg.norm(p1) - L0) / L0 * 100.0)
        out[nm] = float(np.max(dl)) if dl else float("nan")
    return out


print(f"{'pose':>6} {'pen':>5} {'外侧max':>9} {'内侧max':>9} {'外<内?':>8}")
for ie in ("+5", "-3"):
    for pen in ("1", "0.1"):
        o = measure(ie, pen)
        if o is None:
            print(f"{ie:>6} {pen:>5}   构建失败")
            continue
        la = max(o[k] for k in LAT if k in o)
        me = max(o[k] for k in MED if k in o)
        print(f"{ie:>6} {pen:>5} {la:>9.3f} {me:>9.3f} "
              f"{'是 ✅' if la < me else '否 ❌':>8}", flush=True)
