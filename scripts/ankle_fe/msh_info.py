"""解析 gmsh .msh (v2.2) 节点和四面体, 报告数量/体积/边长。

MSH 2.2 element 行格式:
    elm-number  elm-type  number-of-tags  <tag>...  node-number-list
线性四面体 elm-type=4, 节点在 number-of-tags 之后。
"""

import sys
from pathlib import Path

import numpy as np


def parse(path: Path):
    nodes = {}
    tets = []
    section = None
    pending_count = 0
    for raw in path.read_text(errors="replace").splitlines():
        s = raw.strip()
        if s.startswith("$"):
            section = s
            pending_count = 1  # 段头后第一行是计数
            continue
        if not s or pending_count:
            if s and pending_count and section in ("$Nodes", "$Elements"):
                pending_count = 0
                continue
            pending_count = 0 if s else pending_count
            continue
        parts = s.split()
        if section == "$Nodes" and len(parts) >= 4:
            try:
                nodes[int(parts[0])] = (float(parts[1]), float(parts[2]), float(parts[3]))
            except ValueError:
                pass
        elif section == "$Elements" and len(parts) >= 3:
            try:
                etype = int(parts[1])
                ntags = int(parts[2])
            except ValueError:
                continue
            if etype == 4 and len(parts) >= 3 + ntags + 4:
                n = parts[3 + ntags : 3 + ntags + 4]
                tets.append([int(x) for x in n])
    return nodes, tets


for fn in sys.argv[1:]:
    p = Path(fn)
    if not p.exists():
        print(f"[!] 缺 {fn}")
        continue
    nodes, tets = parse(p)
    if not tets:
        print(f"{p.name:<24} 无四面体 (节点 {len(nodes)})")
        continue
    tags = sorted(nodes)
    t2i = {t: i for i, t in enumerate(tags)}
    coord = np.array([nodes[t] for t in tags], dtype=float)
    conn = np.array([[t2i[t] for t in row] for row in tets], dtype=np.int64)
    q = coord[conn]
    v6 = np.einsum(
        "ij,ij->i", q[:, 1] - q[:, 0], np.cross(q[:, 2] - q[:, 0], q[:, 3] - q[:, 0])
    )
    vol = float(np.abs(v6 / 6.0).sum())
    e = np.stack(
        [
            np.linalg.norm(q[:, i] - q[:, j], axis=1)
            for i, j in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
        ],
        axis=1,
    )
    # 单元质量指标: 归一化体积比 (1 = 正四面体)
    qm = 12.0 / np.sqrt(2) * np.abs(v6 / 6.0) / np.power(e.max(1), 3)
    print(
        f"{p.name:<24} 节点 {len(tags):>7,} 四面体 {len(tets):>8,} "
        f"vol={vol*1e6:>8.4f} cm³ | 边 {e.min()*1e3:.3f}~{e.max()*1e3:.2f} mm | "
        f"质量 min/中 {qm.min():.3f}/{np.median(qm):.3f} | 负体积 {(v6<0).sum()}"
    )
