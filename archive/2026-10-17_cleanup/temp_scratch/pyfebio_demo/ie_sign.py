"""独立裁决：POSE_IE=±5 相对 IE=0 的足部组刚体旋转轴/角（判断哪个符号是内翻）。"""
import os
import re
import subprocess
from pathlib import Path

import numpy as np

R = Path(r"D:\Project\climbing_fall_analysis")
W = R / "temp" / "thums" / "febio_lig"
PY = str(R / ".venv" / "Scripts" / "python.exe")
FOOT = ("talus", "calcaneus", "navicular")


def build(ie, out):
    e = dict(os.environ)
    e.update(DRIVER="dispx", DISP_MM="-1.0", LOAD_SIGN="+1",
             FIX_TALUS="transverse", TRUNC_MM="60", AXIS_SRC="full",
             ALIGN_AXIS="1", CART_NEW="1", TIE="1", DET_MIN="0.0",
             TALUS_RIGID="0", TAL_ROTLOCK="0", TIME_STEPS="20", POSE_IE=ie,
             CART_DIR=str(R / "temp" / "thums" / "cart_m1"))
    p = subprocess.run([PY, "-u", str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
                       capture_output=True, text=True, env=e, cwd=str(R), timeout=900)
    so = (p.stdout or "") + (p.stderr or "")
    txt = (W / "ankle_lig_WK.feb").read_text(encoding="latin-1", errors="replace")
    Path(out).write_text(txt, encoding="latin-1")
    # 抓 ML / AP 轴打印
    ml = re.search(r"内外侧轴 ML = \[([^\]]+)\]", so)
    ap = re.search(r"前后轴 AP = \[([^\]]+)\]", so)
    return {int(m.group(1)): [float(m.group(2)), float(m.group(3)),
                              float(m.group(4))]
            for m in re.finditer(r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*'
                                 r'([-\d.eE+]+)\s*,\s*([-\d.eE+]+)', txt)}, so


def nodenames(txt):
    """域 -> 节点 id（从 SolidDomain 段取不到，改从 <NodeSet> 名判断不可靠）
    简易法：用 NodeSet 名含 talus/calcaneus/navicular 的集合。"""
    d = {}
    for m in re.finditer(r'<NodeSet name="([^"]+)">(.*?)</NodeSet>', txt, re.S):
        nm = m.group(1)
        if any(k in nm for k in FOOT):
            d[nm] = [int(x) for x in re.findall(r'<node id="(\d+)"', m.group(2))]
    return d


X0, so0 = build("0", R / "temp" / "pyfebio_demo" / "_ie0.feb")
print("ML =", re.search(r"内外侧轴 ML = \[([^\]]+)\]", so0).group(1))
print("AP =", re.search(r"前后轴 AP = \[([^\]]+)\]", so0).group(1))
txt0 = (R / "temp" / "pyfebio_demo" / "_ie0.feb").read_text(encoding="latin-1",
                                                            errors="replace")
X5, so5 = build("+5", R / "temp" / "pyfebio_demo" / "_ie5.feb")
Xm5, sop = build("-5", R / "temp" / "pyfebio_demo" / "_iem5.feb")

# 用共同节点（节点 id 相同）做 Kabsch
ids = sorted(set(X0) & set(X5))
P0 = np.array([X0[i] for i in ids])
P5 = np.array([X5[i] for i in ids])
Pm = np.array([Xm5[i] for i in ids])
for tag, P1 in (("+5", P5), ("-5", Pm)):
    c0, c1 = P0.mean(0), P1.mean(0)
    H = (P0 - c0).T @ (P1 - c1)
    U, _, Vt = np.linalg.svd(H)
    D = np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))])
    Rm = Vt.T @ D @ U.T
    ang = np.degrees(np.arccos(np.clip((np.trace(Rm) - 1) / 2, -1, 1)))
    w, V = np.linalg.eig(Rm)
    ax = np.real(V[:, np.argmin(np.abs(w - 1))])
    ax = ax / np.linalg.norm(ax)
    print(f"\nPOSE_IE={tag}: 整体刚体旋转 {ang:.3f}°  轴={np.round(ax,4)}")
    print(f"   （ML 轴 = {np.round(np.array([float(x) for x in re.search('ML = .([^]]+).', so0).group(1).split()]),4)}）")

# 只取足部节点（按坐标范围筛：跟骨/距骨/舟骨在关节下方）
zc = np.array([X0[i][2] for i in ids])
sel = zc < np.percentile(zc, 60)
print(f"\n足部子集裁切（按 z 低 60%）：{sel.sum()} / {len(ids)} 节点")
for tag, P1 in (("+5", P5), ("-5", Pm)):
    A, B = P0[sel], P1[sel]
    c0, c1 = A.mean(0), B.mean(0)
    H = (A - c0).T @ (B - c1)
    U, _, Vt = np.linalg.svd(H)
    D = np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))])
    Rm = Vt.T @ D @ U.T
    ang = np.degrees(np.arccos(np.clip((np.trace(Rm) - 1) / 2, -1, 1)))
    w, V = np.linalg.eig(Rm)
    ax = np.real(V[:, np.argmin(np.abs(w - 1))])
    ax = ax / np.linalg.norm(ax)
    print(f"  {tag}: {ang:.3f}°  轴={np.round(ax,4)}")
