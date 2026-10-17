"""读 FEBio log（跳过 *Step 头）。"""
import io
from pathlib import Path

import numpy as np

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio_lig")


def read_cols(fp):
    rows, steps = [], []
    for ln in io.open(fp, encoding="utf-8", errors="replace"):
        ln = ln.strip()
        if ln.startswith("*Step"):
            steps.append(ln)
            continue
        if not ln or ln.startswith("*"):
            continue
        try:
            rows.append([float(x) for x in ln.replace(",", " ").split()])
        except ValueError:
            pass
    return np.array(rows), steps


for nm in ("reac_tibia", "reac_talus", "disp"):
    fs = sorted(W.glob(f"{nm}_*.txt"))
    a, steps = read_cols(fs[0])
    print(f"\n{nm}  {fs[0].name}  {a.shape}  阶段头 {len(steps)}")
    F = np.linalg.norm(a[:, 1:4], axis=1)
    n = len(F)
    print(f"  合力 |F| 曲线（每 5 步）：{np.round(F[::5], 2)}")
    print(f"  最大 {F.max():.2f} @ 行 {F.argmax()}   最小 {F.min():.2f}")
