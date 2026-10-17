"""逐时间点打印支反力三分量：判断 step2 里 z 向是否仍被约束。"""
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


a = rows(next(W.glob("reac_tibia_*.txt")))
nnode = len(np.unique(a[:, 0].astype(int)))
nstep = len(a) // nnode
F = a[:, 1:4].reshape(nstep, nnode, 3).sum(axis=1)
print(f"{nnode} 节点 × {nstep} 时间点；step1 = 前 20 步")
print("  步    Rx        Ry        Rz")
for i, v in enumerate(F):
    tag = "  <-- step1 末" if i == 19 else ("  [step2]" if i >= 20 else "")
    print(f"  {i+1:>3} {v[0]:>9.2f} {v[1]:>9.2f} {v[2]:>9.2f}{tag}")
print("\n施加的载荷（FORCE_N=600，lc=2 曲线在 step2 内 0→1 ⇒ t=1 时 0，"
      "每步 +5%）：")
print("  ", [round(600 * (0.05 * k), 1) for k in range(1, 6)], "N (step2 前 5 步)")
