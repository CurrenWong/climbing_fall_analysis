"""查 prism6 是否"扭转"：底面三角形法向 与 顶面三角形法向 应同向。

扭转(点积<0)的楔形 ⇒ 雅可比为负，但**质心采样可能算不出**（这正是自检漏判的机理）。
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

    tw, flip_bot, badpair = [], [], []
    for k, row in enumerate(E):
        c = [ix(x) for x in row if int(x) != 0]
        if len(c) != 6:
            continue
        T, B = np.array(c[3:]), np.array(c[:3])
        nB = np.cross(N[B[1]] - N[B[0]], N[B[2]] - N[B[0]])
        nT = np.cross(N[T[1]] - N[T[0]], N[T[2]] - N[T[0]])
        cB, cT = N[B].mean(0), N[T].mean(0)
        ax = cT - cB
        dB, dT = float(nB @ ax), float(nT @ ax)
        if dB * dT < 0:
            tw.append(k)
        if dB < 0:
            flip_bot.append(k)
        if not off[B].any() or not off[T].all():
            badpair.append(k)
    print(f"\n===== {j} ({len(E)} prism) =====")
    print(f"  **扭转**(底面法向·顶面法向<0，即 dB*dT<0): {len(tw)} 个  {tw[:12]}")
    print(f"  底面法向与轴反向 (dB<0): {len(flip_bot)} 个  {flip_bot[:12]}")
    print(f"  底面/顶面 offset 约定不符: {len(badpair)} 个")
