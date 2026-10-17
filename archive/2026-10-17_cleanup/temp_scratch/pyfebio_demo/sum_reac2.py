"""正确重塑：165 节点 × 25 时间点；每时间点对节点求和。"""
import io
from pathlib import Path

import numpy as np

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio_lig")


def rows(fp):
    out = []
    for ln in io.open(fp, encoding="utf-8", errors="replace"):
        ln = ln.strip()
        if not ln or ln.startswith("*"):
            continue
        try:
            out.append([float(x) for x in ln.replace(",", " ").split()])
        except ValueError:
            pass
    return np.array(out)


for nm in ("reac_tibia", "reac_fibula", "reac_talus"):
    a = rows(next(W.glob(f"{nm}_*.txt")))
    nid = a[:, 0].astype(int)
    nnode = len(np.unique(nid))
    nstep = len(a) // nnode
    F = a[:, 1:4].reshape(nstep, nnode, 3).sum(axis=1)
    mag = np.linalg.norm(F, axis=1)
    print(f"\n{nm}: {nnode} 节点 × {nstep} 时间点")
    print("  逐步 合力 N（20=step1末）：")
    print("   ", np.round(mag, 2))
    if nstep >= 20:
        print(f"  ⇒ 第 20 步(step1 末) 合力 = {mag[19]:.2f} N，向量 {np.round(F[19],2)}")
