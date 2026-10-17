"""反力落在哪些节点？→ 揭示真实载荷路径。"""
import re
from pathlib import Path

import h5py
import numpy as np

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio")
AX = np.array([0.8153, 0.135, 0.5631])
AX = AX / np.linalg.norm(AX)

for tag, h5 in (("B(基线 PENALTY=1)", "_pB.h5"), ("C(基线 AUGLAG)", "_pC.h5")):
    txt = (W / f"_{tag[0]}.feb").read_text(encoding="latin-1")
    parts = {}
    for m in re.finditer(
            r'<Elements\b[^>]*type="(tet4|hex8|penta6)"[^>]*name="([^"]+)"[^>]*>(.*?)</Elements>',
            txt, re.S):
        parts[m.group(2)] = {int(y) for s in
                             re.findall(r"<elem\b[^>]*>\s*([\d,]+)\s*</elem>", m.group(3))
                             for y in s.split(",")}
    nss = {m.group(1): {int(x) for x in m.group(2).split(",") if x.strip()}
           for m in re.finditer(r'<NodeSet name="([^"]+)">([^<]*)</NodeSet>', txt)}
    bcs = re.findall(r'<bc\b[^>]*name="([^"]+)"[^>]*node_set="([^"]+)"[^>]*>', txt)
    with h5py.File(W / h5, "r") as h:
        ns = sorted((int(q) for q in h["states"]), key=int)
        R = np.asarray(h[f"states/{ns[-1]}/node_data/reaction forces/all"][()], float)
        U = np.asarray(h[f"states/{ns[-1]}/node_data/displacement/all"][()], float)
    rn = np.linalg.norm(R, axis=1)
    nz = np.where(rn > 1e-6)[0] + 1
    print("=" * 78)
    print(f"{tag}   非零反力节点 {len(nz)}/{len(R)}")
    print("=" * 78)
    print(f"  BC: {bcs}")
    for k, v in nss.items():
        print(f"  NodeSet {k:26s} {len(v)}")
    tot = (R @ AX).sum()
    print(f"  沿轴支反力总和 = {tot:.2f} N")
    # 非零反力节点的归属
    lab = {}
    for i in nz:
        got = None
        for d, s in parts.items():
            if i in s:
                got = d
        lab.setdefault(got, []).append(i)
    print("  非零反力节点的部件归属：")
    for d, v in sorted(lab.items(), key=lambda kv: -len(kv[1])):
        v = np.array(v)
        rr = R[v - 1] @ AX
        print(f"    {str(d):28s} {len(v):>5d} 节点  沿轴反力 {rr.sum():>9.2f} N  "
              f"|R| 中位 {np.median(np.linalg.norm(R[v-1],axis=1)):.4f}")
    # 与 BC 集合的交叉
    for k, v in nss.items():
        inter = set(nz) & v
        if inter:
            vv = np.array(sorted(inter))
            print(f"    ∩ {k:24s} {len(inter):>5d}  沿轴反力 "
                  f"{(R[vv-1] @ AX).sum():>9.2f} N")
    # 非BC节点为何有反力？
    allbc = set().union(*nss.values()) if nss else set()
    extra = set(nz) - allbc
    print(f"  ⚠ 不在任何 NodeSet 里却有反力的节点 = {len(extra)}")
    if extra:
        ev = np.array(sorted(extra))
        print(f"     它们的沿轴反力 = {(R[ev-1] @ AX).sum():.4f} N "
              f"|R| 中位 {np.median(np.linalg.norm(R[ev-1],axis=1)):.6f} "
              f"max {np.linalg.norm(R[ev-1],axis=1).max():.4f}")
