"""按时间步分组求和：step1 末尾的总支反力到底多大。"""
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


for nm in ("reac_tibia", "reac_talus"):
    a = rows(next(W.glob(f"{nm}_*.txt")))
    t = a[:, 0]
    # 按 time 值分组
    tot = []
    for tv in np.unique(t):
        blk = a[np.isclose(t, tv)]
        tot.append((tv, np.linalg.norm(blk[:, 1:4].sum(axis=0))))
    tot = np.array(tot)
    print(f"\n{nm}: {len(tot)} 个时间点（每点 {int(np.sum(np.isclose(t, t[0])))} 节点）")
    print("  t / 合力 N:")
    for tv, F in tot[::max(1, len(tot) // 22)]:
        tag = "  <-- step1 末" if abs(tv - 1.0) < 1e-6 else ""
        print(f"    {tv:8.5f}  {F:10.3f}{tag}")
    i1 = int(np.argmin(np.abs(tot[:, 0] - 1.0)))
    print(f"  ⇒ t=1 处合力 = {tot[i1,1]:.3f} N ；t 最大 {tot[-1,0]:.4f} "
          f"合力 {tot[-1,1]:.3f} N")
