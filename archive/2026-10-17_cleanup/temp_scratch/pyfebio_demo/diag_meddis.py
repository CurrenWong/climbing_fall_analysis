"""我的 cart_patch(m=1) 几何 vs 基线：外层节点位置差多少？
若同一骨节点的各面 offset 差异大 ⇒ 用「逐节点中位数」会把几何改掉。
"""
import numpy as np
from pathlib import Path

THUMS = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
src = np.load(THUMS / "cartilage_mesh.npz", allow_pickle=True)
bpos = {int(v): src["bone_xyz"][i] for i, v in enumerate(src["bone_node_ids"])}

for j in ["tibiotalar_tibia", "tibiotalar_talus",
          "subtalar_talus", "subtalar_calcaneus"]:
    E = np.asarray(src[f"{j}_elems"], np.int64)
    off = src[f"{j}_offsets"]
    dacc: dict[int, list[np.ndarray]] = {}
    pairs = []
    for row in E:
        v = [int(x) for x in row if int(x) != 0]
        m = len(v) // 2
        for r, o in zip(v[:m], v[m:]):
            d = off[-o - 1] - bpos[r]
            dacc.setdefault(r, []).append(d)
            pairs.append((r, -o - 1))

    nbase = len(dacc)
    nouter = len({o for _, o in pairs})
    # 每个骨节点被几个面用、每次用的 outer 是否同一个
    n_outer_per_base = {r: len({o for rr, o in pairs if rr == r}) for r in dacc}
    multi = sum(1 for r, k in n_outer_per_base.items() if k > 1)
    # 中位数 vs 各自实际位移的偏差
    dev = []
    for r, ds in dacc.items():
        med = np.median(np.array(ds), axis=0)
        for d in ds:
            dev.append(np.linalg.norm(d - med))
    dev = np.array(dev)
    # 同一节点内各方向的夹角（>30° 说明方向不一致）
    angs = []
    for r, ds in dacc.items():
        if len(ds) < 2:
            continue
        A = np.array(ds)
        for i in range(len(A)):
            for k2 in range(i + 1, len(A)):
                c = np.dot(A[i], A[k2]) / (np.linalg.norm(A[i]) * np.linalg.norm(A[k2]) + 1e-15)
                angs.append(np.degrees(np.arccos(np.clip(c, -1, 1))))
    angs = np.array(angs) if angs else np.array([0.0])
    print(f"[{j:22s}] 底节点 {nbase:4d}  外节点 {nouter:4d}  "
          f"一底多外 {multi:3d}  中位数vs实际 偏差中位 {np.median(dev):.4f} mm "
          f"最大 {dev.max():.4f}  同节点方向夹角 中位 {np.median(angs):5.1f}° "
          f"最大 {angs.max():5.1f}°")
