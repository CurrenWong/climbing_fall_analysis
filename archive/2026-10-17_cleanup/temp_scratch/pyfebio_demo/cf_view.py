"""CARTFIX 诊断 8：把底面补丁画出来（含坏单元高亮 + disp 方向箭头）。"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
CM1 = np.load(ROOT / "temp/thums/cart_m1/cartilage_mesh.npz", allow_pickle=True)
RS = (0.0, 1 / 6, 1 / 3, 1 / 2, 2 / 3, 1.0)
TS = (0.0, 0.25, 0.5, 0.75, 1.0)
GS = [(r, s, t) for r in RS for s in RS if r + s <= 1 + 1e-9 for t in TS]


def dets(P6):
    out = np.full(len(P6), np.inf)
    for (r, s, t) in GS:
        G = np.array([
            [-(1 - t), -(1 - t), -(1 - r - s)],
            [1 - t, 0.0, -r],
            [0.0, 1 - t, -s],
            [-t, -t, (1 - r - s)],
            [t, 0.0, r],
            [0.0, t, s],
        ])
        out = np.minimum(out, np.linalg.det(np.einsum("nia,ic->nac", P6, G)))
    return out


def ix(x):
    x = int(x)
    return (-x - 1) if x < 0 else (x - 1)


fig = plt.figure(figsize=(16, 10))
for col, j in enumerate(["tibiotalar_tibia", "tibiotalar_talus"]):
    N = np.asarray(CM1[f"{j}_nodes"], float)
    E = np.asarray(CM1[f"{j}_elems"], np.int64)
    C = np.array([[ix(x) for x in row] for row in E])
    B, T = C[:, :3], C[:, 3:]
    W = dets(N[C])
    d = N[T] - N[B]
    ax = fig.add_subplot(2, 2, 1 + col * 2, projection="3d")
    for i in range(len(E)):
        P = N[B[i]]
        bad = W[i] <= 1e-2
        tri = np.vstack([P, P[:1]])
        ax.plot(tri[:, 0], tri[:, 1], tri[:, 2],
                color=("red" if bad else "steelblue"), lw=(1.6 if bad else 0.6),
                alpha=1.0 if bad else 0.5)
    allb = N[B.ravel()]
    ctr = allb.mean(0)
    ax.set_title(f"{j} 底面补丁 (红=detJ≤0.01)  n={len(E)}")
    ax.set_box_aspect((1, 1, 1))
    # 第二个子图：disp 方向 + 顶面
    ax2 = fig.add_subplot(2, 2, 2 + col * 2, projection="3d")
    for i in range(0, len(E), 1):
        P = N[B[i]]
        bad = W[i] <= 1e-2
        tri = np.vstack([P, P[:1]])
        ax2.plot(tri[:, 0], tri[:, 1], tri[:, 2],
                 color=("red" if bad else "steelblue"), lw=(1.4 if bad else 0.5),
                 alpha=1.0 if bad else 0.35)
    for i in range(len(E)):
        for a in range(3):
            b, t = int(B[i, a]), int(T[i, a])
            ax2.plot([N[b, 0], N[t, 0]], [N[b, 1], N[t, 1]], [N[b, 2], N[t, 2]],
                     color="orange", lw=0.5, alpha=0.7)
    ax2.set_title(f"{j} 底面+偏置向量(橙)  中心 {np.round(ctr,1)}")
    ax2.set_box_aspect((1, 1, 1))
    print(f"{j}: 补丁质心 {np.round(ctr,2)}  bbox {np.round(allb.min(0),1)} .. {np.round(allb.max(0),1)}")
fig.tight_layout()
out = ROOT / "temp/pyfebio_demo/cf_patch_view.png"
fig.savefig(out, dpi=95)
print("saved", out)
