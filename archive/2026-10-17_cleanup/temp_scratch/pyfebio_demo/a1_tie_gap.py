"""测 tie 的 primary(软骨底面节点) 到 secondary(骨面) 的实际间隙 —— 两侧对比。
（正对照 = tibiotalar 两对：已在跑、能收敛）"""
import re
from pathlib import Path

import numpy as np

feb = Path(r"D:\Project\climbing_fall_analysis\temp/thums/febio_lig/ankle_lig_WK.feb")
t = feb.read_text(encoding="latin-1")
xyz = {}
for m in re.finditer(r'<node id="(\d+)">([-0-9.eE+,]+)</node>', t):
    xyz[int(m.group(1))] = np.array([float(x) for x in m.group(2).split(",")])


def surf(name):
    m = re.search(rf'<Surface name="{name}">(.*?)</Surface>', t, re.S)
    if not m:
        return []
    return [[int(x) for x in g.split(",")]
            for g in re.findall(r'<(?:tri3|quad4) id="\d+">([0-9,]+)</(?:tri3|quad4)>',
                                m.group(1))]


def pt_tri(p, T):
    a, b, c = T
    ab, ac = b - a, c - a
    ap = p - a
    d1, d2 = ab @ ap, ac @ ap
    bp = p - b
    d3, d4 = ab @ bp, ac @ bp
    cp = p - c
    d5, d6 = ab @ cp, ac @ cp
    vc = d1 * d4 - d3 * d2
    vb = d5 * d2 - d1 * d6
    va = d3 * d6 - d5 * d4
    if d1 <= 0 and d2 <= 0:
        return np.linalg.norm(ap)
    if d3 >= 0 and d4 <= d3:
        return np.linalg.norm(bp)
    if vc <= 0 and d1 >= 0 and d3 <= 0:
        v = d1 / (d1 - d3)
        return np.linalg.norm(a + v * ab - p)
    if d6 >= 0 and d5 <= d6:
        return np.linalg.norm(cp)
    if vb <= 0 and d2 >= 0 and d6 <= 0:
        w = d2 / (d2 - d6)
        return np.linalg.norm(a + w * ac - p)
    if va <= 0 and (d4 - d3) >= 0 and (d5 - d6) >= 0:
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        return np.linalg.norm(b + w * (c - b) - p)
    den = 1.0 / (va + vb + vc)
    v, w = vb * den, vc * den
    return np.linalg.norm(a + ab * v + ac * w - p)


def report(prim_name, sec_name):
    P = surf(prim_name)
    S = surf(sec_name)
    if not P or not S:
        print(f"{prim_name} / {sec_name}: 缺失")
        return
    nodes = sorted({x for f in P for x in f})
    tris = []
    for f in S:
        if len(f) == 3:
            tris.append([xyz[i] for i in f])
        else:
            tris.append([xyz[f[0]], xyz[f[1]], xyz[f[2]]])
            tris.append([xyz[f[0]], xyz[f[2]], xyz[f[3]]])
    ST = np.array(tris)
    # 预算：用 secondary 面片质心做粗筛，再精确
    cen = ST.mean(1)
    d = []
    for n in nodes:
        p = xyz[n]
        k = int(np.argmin(np.linalg.norm(cen - p, axis=1)))
        # 粗筛后在 6 个近邻面片里精算
        order = np.argsort(np.linalg.norm(cen - p, axis=1))[:8]
        d.append(min(pt_tri(p, ST[j]) for j in order))
    d = np.array(d)
    print(f"{prim_name} ({len(nodes)} 节点) → {sec_name} ({len(S)} 面):")
    print(f"   间隙 min/p5/median/p95/max = {d.min():.4f}/{np.percentile(d,5):.4f}"
          f"/{np.median(d):.4f}/{np.percentile(d,95):.4f}/{d.max():.4f} mm"
          f"   >0.01mm: {int((d>0.01).sum())}  >0.05mm: {int((d>0.05).sum())}"
          f"  >0.5mm: {int((d>0.5).sum())}")


report("tiesurf_tibiotalar_tibia", "bonesurf_tibiotalar_tibia")
report("tiesurf_tibiotalar_talus", "bonesurf_tibiotalar_talus")
report("tiesurf_subtalar_calcaneus", "bonesurf_subtalar_calcaneus")
