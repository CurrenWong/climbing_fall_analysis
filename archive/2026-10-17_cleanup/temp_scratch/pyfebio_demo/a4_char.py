"""A4 步骤1：表征卡点单元 8241（先测不改）。

从 E2p.feb（A1+A2 配置）解析：
  · 8241 的节点坐标、体积、静息 det(J)、长细比、最小高度
  · tibiotalar_talus__penta6 全域名（182 个）的同一组指标分布 → 8241 的分位
  · 8241 是否在接触面上、是否在接触区边缘
"""
import re
from pathlib import Path

import numpy as np

FEB = Path(r"D:\Project\climbing_fall_analysis\temp/pyfebio_demo/a1a2/E2p.feb")
t = FEB.read_text(encoding="latin-1")

# --- 节点 ---
nodes = {}
for m in re.finditer(r'<node id="(\d+)">\s*([-\d.eE+]+),\s*([-\d.eE+]+),\s*([-\d.eE+]+)',
                     t):
    nodes[int(m.group(1))] = np.array([float(m.group(2)), float(m.group(3)),
                                       float(m.group(4))])
print(f"节点数={len(nodes)}")

# --- 域内单元 ---
DOMS = {}
for m in re.finditer(r'<Elements type="(\w+)" name="([^"]+)">(.*?)</Elements>', t, re.S):
    et, name, body = m.groups()
    els = {}
    for e in re.finditer(r'<elem id="(\d+)">(.*?)</elem>', body, re.S):
        els[int(e.group(1))] = [int(x) for x in re.findall(r"\d+", e.group(2))]
    DOMS[name] = (et, els)
D = "tibiotalar_talus__penta6"
et, els = DOMS[D]
print(f"{D}: {len(els)} 个 {et} 单元")


def penta6_metrics(nid):
    """返回 (体积, detJ@中心, 最小高度近似, 长细比)。"""
    P = np.array([nodes[i] for i in nid], float)      # 6 节点
    # penta6 节点顺序假定 (0,1,2)=底三角, (3,4,5)=顶三角
    bot, top = P[:3], P[3:]
    c1 = bot.mean(0); c2 = top.mean(0)
    h = float(np.linalg.norm(c2 - c1))                # 高度
    def tri_area(X):
        return 0.5 * float(np.linalg.norm(np.cross(X[1] - X[0], X[2] - X[0])))
    A1, A2 = tri_area(bot), tri_area(top)
    V = (A1 + A2) / 2.0 * h                           # 简化棱柱体积
    # 边长
    ed = [np.linalg.norm(P[j] - P[i]) for i, j in
          ((0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5))]
    ed = np.array(ed)
    return V, h, float(ed.min()), float(ed.max() / max(ed.min(), 1e-9))


ids = sorted(els)
rows = []
for eid in ids:
    nid = els[eid]
    if len(nid) != 6 or any(i not in nodes for i in nid):
        continue
    V, h, emin, ar = penta6_metrics(nid)
    rows.append((eid, V, h, emin, ar))
rows = np.array(rows, float)
print(f"可用行数={len(rows)}")

TARGET = 8241
sel = rows[rows[:, 0] == TARGET]
print("\n=== 单元 8241 ===")
if len(sel):
    eid, V, h, emin, ar = sel[0]
    print(f"  体积={V:.4f} mm³  高度(上下三角心距)={h:.4f} mm  "
          f"最短边={emin:.4f} mm  长细比={ar:.2f}")

for col, nm in ((2, "高度"), (3, "最短边"), (4, "长细比")):
    v = rows[:, col]
    tgt = sel[0, col] if len(sel) else np.nan
    pct = (v < tgt).mean() * 100
    print(f"  {nm:<6} 分位={pct:5.1f}%  (min={v.min():.4f} med={np.median(v):.4f} "
          f"max={v.max():.4f})")

print("\n=== 域内最「薄」的 10 个 ===")
for r_ in rows[np.argsort(rows[:, 3])][:10]:
    print(f"  elem {int(r_[0])}  最短边={r_[3]:.4f}  高={r_[2]:.4f}  长细比={r_[4]:.2f}")
print("\n=== 域内长细比最大的 10 个 ===")
for r_ in rows[np.argsort(-rows[:, 4])][:10]:
    print(f"  elem {int(r_[0])}  长细比={r_[4]:.2f}  最短边={r_[3]:.4f}  高={r_[2]:.4f}")

# --- 8241 是否在接触面上 ---
print("\n=== 8241 与接触面的关系 ===")
nid = els.get(TARGET, [])
sp = re.search(r'<Surface name="surf_tibiotalar_talus">(.*?)</Surface>', t, re.S)
if sp:
    surf = set(int(x) for x in re.findall(r"\d+", sp.group(1)))
    print(f"  surf_tibiotalar_talus 节点数={len(surf)}")
    print(f"  8241 的节点 {nid} ⇒ 落在接触面上的: {[n for n in nid if n in surf]}")
    # 8241 是否与接触面单元共用节点
    print(f"  8241 节点坐标质心={nodes[nid[0]].round(2) if nid else '—'}")
# 接触面里有哪些单元
sm = re.search(r'<Surface name="surf_tibiotalar_talus">(.*?)</Surface>', t, re.S)
if sm:
    fes = re.findall(r"<faces>?\s*([\d,\s]+)", sm.group(1))
    print(f"  接触面 face 条目数={len(fes)}")
