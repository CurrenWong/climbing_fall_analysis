"""查 ATFL（前距腓）为何在内翻下几乎不受力：看它的节点对与坐标。"""
import io
import os
import re
import subprocess
from pathlib import Path

import numpy as np

R = Path(r"D:\Project\climbing_fall_analysis")
W = R / "temp" / "thums" / "febio_lig"
PY = str(R / ".venv" / "Scripts" / "python.exe")
e = dict(os.environ)
e.update(DRIVER="dispx", DISP_MM="-1.0", LOAD_SIGN="+1", FIX_TALUS="transverse",
         TRUNC_MM="60", AXIS_SRC="full", ALIGN_AXIS="1", CART_NEW="1", TIE="1",
         DET_MIN="0.0", TALUS_RIGID="0", TAL_ROTLOCK="0", TIME_STEPS="20",
         POSE_IE="+5", LOG_DATA="1",
         CART_DIR=str(R / "temp" / "thums" / "cart_m1"))
subprocess.run([PY, "-u", str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
               capture_output=True, text=True, env=e, cwd=str(R), timeout=900)
txt = (W / "ankle_lig_WK.feb").read_text(encoding="latin-1", errors="replace")
X = {int(m.group(1)): np.array([float(m.group(2)), float(m.group(3)),
                                float(m.group(4))])
     for m in re.finditer(r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)'
                          r'\s*,\s*([-\d.eE+]+)', txt)}
z = np.load(W / "lig_read_WK.npz", allow_pickle=True)
names = [str(x) for x in z["lig_names"]]
PR, npr, ea = z["lig_pairs"], z["lig_npair"], z["lig_ea"]
# 韧带的 FE 节点归属（韧带块的单元 → 节点）
lig_elem = {}
for m in re.finditer(r'<Elements type="([a-z0-9]+)" name="([^"]+)">(.*?)'
                     r'</Elements>', txt, re.S):
    nm = m.group(2)
    ids = [int(y) for x in re.findall(r'<elem id="\d+">\s*([\d,]+)\s*</elem>',
                                      m.group(3)) for y in x.split(",")]
    if ids:
        lig_elem[nm] = ids
print("韧带元素块:", {k: len(v) for k, v in list(lig_elem.items())[:14]})
for i, nm in enumerate(names):
    k = int(npr[i])
    Ls = []
    for a, b in PR[i, :k, :]:
        a, b = int(a), int(b)
        if a < 0 or b < 0 or a not in X or b not in X:
            continue
        Ls.append(float(np.linalg.norm(X[a] - X[b])))
    if not Ls:
        continue
    print(f"\n{nm}  EA={ea[i]:.1f}  n={len(Ls)}")
    print(f"   长度 mm: min {min(Ls):.2f}  中位 {np.median(Ls):.2f}  max {max(Ls):.2f}")
    pr0 = list(PR[i, :min(k, 3), :])
    for a, b in pr0:
        a, b = int(a), int(b)
        if a in X and b in X:
            print(f"   pair ({a},{b})  L={np.linalg.norm(X[a]-X[b]):.2f} mm  "
                  f"Δz={X[a][2]-X[b][2]:+.2f}")
