"""量 cart_patch 输出里每个 prism6 的厚度（底面质心 -> 顶面质心）与扁的程度。

目的：找出让 FEBio 报「2 negative jacobians」的那 2 个极扁单元。
用法: python prism_thick.py
"""
from pathlib import Path

import numpy as np

P = Path(r"D:\Project\climbing_fall_analysis\temp\thums\cart_new\cartilage_mesh.npz")
d = np.load(P)
for j in ("tibiotalar_tibia", "tibiotalar_talus"):
    N = np.asarray(d[f"{j}_nodes"], float)
    E = np.asarray(d[f"{j}_elems"])
    off = np.asarray(d[f"{j}_is_offset"], bool)

    def ix(x):
        x = int(x)
        return (-x - 1) if x < 0 else (x - 1)

    print(f"\n===== {j}  ({len(E)} prism, {len(N)} nodes) =====")
    bad = 0
    rows = []
    for k, row in enumerate(E):
        c = [ix(x) for x in row if int(x) != 0]
        if len(c) != 6:
            rows.append((k, np.nan, np.nan, 0))
            continue
        bot, top = np.array(c[:3]), np.array(c[3:])
        cb, ct = N[bot].mean(0), N[top].mean(0)
        h = float(np.linalg.norm(ct - cb))
        # 面积（底面三角形）
        A = 0.5 * float(np.linalg.norm(np.cross(N[bot[1]] - N[bot[0]], N[bot[2]] - N[bot[0]])))
        # 底面 3 节点是否全为非 offset；顶面是否全 offset
        okb = not off[bot].any()
        okt = off[top].all()
        rows.append((k, h, A, 1 if (okb and okt) else 0))
        if not (okb and okt):
            bad += 1
    hs = np.array([r[1] for r in rows if r[1] == r[1]])
    As = np.array([r[2] for r in rows if r[2] == r[2]])
    print(f"  底面/顶面 offset 约定不符的 prism: {bad}")
    print(f"  厚度 h: min={hs.min():.5f}  p1={np.percentile(hs,1):.5f}  中位={np.median(hs):.5f}  "
          f"p99={np.percentile(hs,99):.5f}  max={hs.max():.5f} mm")
    print(f"  底面积 A: min={As.min():.6f}  中位={np.median(As):.6f} mm^2")
    order = np.argsort(hs)
    print("  最扁 8 个（按厚度）:")
    for t in order[:8]:
        k, h, A, ok = rows[int(t)]
        print(f"    prism#{k:<5} h={h:.6f} mm  A={A:.6f} mm^2  h/A^0.5={h/max(A,1e-12)**0.5:.4f}  约定={'OK' if ok else 'BAD'}")
    thin = int((hs < 0.05).sum())
    print(f"  厚度 < 0.05 mm 的 prism: {thin} 个")
    print(f"  厚度 < 0.02 mm 的 prism: {int((hs < 0.02).sum())} 个")
