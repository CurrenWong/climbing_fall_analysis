"""调试：cartilage_mesh.npz 的真实结构（单元索引基准、节点含义）。"""
from pathlib import Path

import numpy as np

R = Path(r"D:\Project\climbing_fall_analysis")
for tag in ("cart_m1",):
    z = np.load(R / "temp/thums" / tag / "cartilage_mesh.npz", allow_pickle=True)
    print("=" * 96)
    print(f"{tag}  键 {len(z.files)}")
    print("=" * 96)
    for k in sorted(z.files):
        v = z[k]
        try:
            print(f"  {k:<34} shape={str(v.shape):<16} dtype={v.dtype}")
        except Exception:
            print(f"  {k:<34} {type(v)}")
    j = "tibiotalar_talus"
    N = np.asarray(z[f"{j}_nodes"], float)
    E = np.asarray(z[f"{j}_elems"], int)
    print(f"\n--- {j} ---")
    print(f"  nodes {N.shape}  x[{N[:,0].min():.1f},{N[:,0].max():.1f}] "
          f"y[{N[:,1].min():.1f},{N[:,1].max():.1f}] "
          f"z[{N[:,2].min():.1f},{N[:,2].max():.1f}]")
    print(f"  elems {E.shape}  min={E.min()} max={E.max()}")
    print(f"  ⇒ 索引基准：min={E.min()} ⇒ {'1-based' if E.min()>=1 else '0-based'}")
    print(f"  前 3 行元素：\n{E[:3]}")
    # 试两种基准各算一个单元的尺寸
    for base in (0, 1):
        try:
            e = E[0] - base
            if e.min() < 0 or e.max() >= len(N):
                print(f"  base={base}: 越界，不可行")
                continue
            P = N[e]
            bot = P[[0, 1, 2]].mean(0)
            top = P[[3, 4, 5]].mean(0)
            print(f"  base={base}: 厚度 {np.linalg.norm(top-bot):.4f} mm  "
                  f"底三角边 "
                  f"{[round(float(np.linalg.norm(P[i]-P[j])),3) for i,j in ((0,1),(1,2),(2,0))]}"
                  f"  顶三角边 "
                  f"{[round(float(np.linalg.norm(P[i]-P[j])),3) for i,j in ((3,4),(4,5),(5,3))]}")
        except Exception as ex:
            print(f"  base={base}: {ex}")
    # 其他键看看
    for k in sorted(z.files):
        if k.endswith("_is_offset") or k.endswith("_tie_faces"):
            v = z[k]
            print(f"  [{k}] shape {getattr(v,'shape','?')} 前几个 {np.asarray(v).ravel()[:8]}")
