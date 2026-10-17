"""校验 _tri_closest 的「内部」判据对三角绕向是否不变（镜像三角形）。"""
import io
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
SRC = io.open(ROOT / "scripts/ankle_fe/cart_patch.py", encoding="utf-8").read()
ns = {"__name__": "h"}
exec(compile(SRC.split("\nsrc = np.load(")[0], "head", "exec"), ns)
tc = ns["_tri_closest"]

rng = np.random.default_rng(0)
T = rng.normal(size=(300, 3, 3))
P = rng.normal(size=(300, 3))
# 镜像：交换 b,c
Tm = T[:, [0, 2, 1]]
q1, i1 = tc(P, T[:, 0], T[:, 1], T[:, 2])
q2, i2 = tc(P, Tm[:, 0], Tm[:, 1], Tm[:, 2])
d1 = np.linalg.norm(q1 - P, axis=1)
d2 = np.linalg.norm(q2 - P, axis=1)
print("距离一致?", np.allclose(d1, d2, atol=1e-9), " 最大差", np.abs(d1 - d2).max())
print("内部标志一致?", bool((i1 == i2).all()), " 不一致数", int((i1 != i2).sum()))
# 与「暴力采样」对比最近距离
S = np.concatenate([T[:, i] for i in range(3)] + [(T[:, 0] + T[:, 1]) / 2,
                                                 (T[:, 1] + T[:, 2]) / 2,
                                                 (T[:, 2] + T[:, 0]) / 2,
                                                 T.mean(1)], axis=1)
print(f"（粗对照）解析最近距 ≤ 采样最近距? {bool((d1 <= np.array([np.linalg.norm(P[i][None]-S[i],axis=1).min() for i in range(len(P))])+1e-9).all())}")
