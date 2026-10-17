"""可行性：关节底面补丁能否用「最佳拟合平面」参数化？
- 浅度 = 平面外最大偏离 / 补丁直径（越小越好；>0.5 说明弯得很厉害）
- 折叠检查 = 投影后每个面的**有向面积符号**是否翻转（翻转 ⇒ 投影不可逆）
- 边界环数 = 1 才算磁盘状（可参数化）
"""
import numpy as np
from pathlib import Path
from collections import Counter, defaultdict

THUMS = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
src = np.load(THUMS / "cartilage_mesh.npz", allow_pickle=True)
bpos = {int(v): src["bone_xyz"][i] for i, v in enumerate(src["bone_node_ids"])}

print("=" * 92)
print(f"{'关节':22s} {'节点':>5s} {'面':>5s} {'tri':>4s} {'边界环':>6s} "
      f"{'边界边':>6s} {'浅度':>7s} {'面积比':>7s} {'翻转面':>7s}")
print("=" * 92)
for j in ["tibiotalar_tibia", "tibiotalar_talus",
          "subtalar_talus", "subtalar_calcaneus"]:
    E = np.asarray(src[f"{j}_elems"], np.int64)
    faces = []
    for row in E:
        v = [int(x) for x in row if int(x) != 0]
        m = len(v) // 2
        faces.append(v[:m])

    nodes = sorted({r for f in faces for r in f})
    P = np.array([bpos[r] for r in nodes])
    c = P.mean(axis=0)
    # PCA 平面
    _, _, Vt = np.linalg.svd(P - c, full_matrices=False)
    n = Vt[2]
    e1, e2 = Vt[0], Vt[1]

    dev = np.abs((P - c) @ n)
    diam = 2 * np.sqrt(np.max(((P - c) ** 2).sum(axis=1)))
    shallow = dev.max() / diam

    # 投影
    UV = np.stack([(P - c) @ e1, (P - c) @ e2], axis=1)
    uv = {r: UV[i] for i, r in enumerate(nodes)}

    def sarea(p2):
        s = 0.0
        for i in range(len(p2)):
            a, b = p2[i], p2[(i + 1) % len(p2)]
            s += a[0] * b[1] - b[0] * a[1]
        return 0.5 * s

    nflip = 0
    a3 = a2 = 0.0
    for f in faces:
        p3 = np.array([bpos[r] for r in f])
        v3 = np.zeros(3)
        for i in range(len(p3)):
            v3 += np.cross(p3[i], p3[(i + 1) % len(p3)])
        a3 += 0.5 * np.linalg.norm(v3)
        a2 += abs(sarea(np.array([uv[r] for r in f])))
        if sarea(np.array([uv[r] for r in f])) < 0:
            nflip += 1

    # 边界环
    ec = Counter()
    for f in faces:
        for i in range(len(f)):
            a, b = f[i], f[(i + 1) % len(f)]
            ec[(min(a, b), max(a, b))] += 1
    bnd = [k for k, v in ec.items() if v == 1]
    # 连通分量数
    adj = defaultdict(set)
    for a, b in bnd:
        adj[a].add(b)
        adj[b].add(a)
    seen, loops = set(), 0
    for a in list(adj):
        if a in seen:
            continue
        loops += 1
        st = [a]
        while st:
            x = st.pop()
            if x in seen:
                continue
            seen.add(x)
            st.extend(adj[x] - seen)

    print(f"{j:22s} {len(nodes):5d} {len(faces):5d} "
          f"{sum(1 for f in faces if len(f)==3):4d} {loops:6d} {len(bnd):6d} "
          f"{shallow:7.3f} {a2/max(a3,1e-9):7.3f} {nflip:7d}")
print("=" * 92)
print("判据：浅度 <0.3 且 翻转面=0 且 边界环=1 ⇒ 平面参数化可用")
