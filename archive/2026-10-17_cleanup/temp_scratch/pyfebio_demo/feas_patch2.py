"""修正版可行性：① 强制 (e1,e2,n) 右手 ② 连通分量 ③ 真折叠数 ④ 边界面构成"""
import numpy as np
from pathlib import Path
from collections import Counter, defaultdict

THUMS = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
src = np.load(THUMS / "cartilage_mesh.npz", allow_pickle=True)
bpos = {int(v): src["bone_xyz"][i] for i, v in enumerate(src["bone_node_ids"])}

print("=" * 96)
print(f"{'关节':22s} {'节点':>5s} {'面':>5s} {'连通块':>6s} {'块大小':>18s} "
      f"{'浅度':>7s} {'真折叠':>6s} {'边界环':>6s} {'边界占比':>8s}")
print("=" * 96)
for j in ["tibiotalar_tibia", "tibiotalar_talus",
          "subtalar_talus", "subtalar_calcaneus"]:
    E = np.asarray(src[f"{j}_elems"], np.int64)
    faces = []
    for row in E:
        v = [int(x) for x in row if int(x) != 0]
        faces.append(v[: len(v) // 2])

    nodes = sorted({r for f in faces for r in f})
    P = np.array([bpos[r] for r in nodes])
    c = P.mean(axis=0)
    _, _, Vt = np.linalg.svd(P - c, full_matrices=False)
    n, e1, e2 = Vt[2], Vt[0], Vt[1]
    # ★ 强制右手系：否则所有面投影后符号整体翻转（假折叠）
    if np.dot(np.cross(e1, e2), n) < 0:
        e2 = -e2
    UV = np.stack([(P - c) @ e1, (P - c) @ e2], axis=1)
    uv = {r: UV[i] for i, r in enumerate(nodes)}

    def sarea(p2):
        s = 0.0
        for i in range(len(p2)):
            a, b = p2[i], p2[(i + 1) % len(p2)]
            s += a[0] * b[1] - b[0] * a[1]
        return 0.5 * s

    # 面法向 vs 平面法向的符号（判断该面是否"朝上"）
    def fnormal(p3):
        v = np.zeros(3)
        for i in range(len(p3)):
            v += np.cross(p3[i], p3[(i + 1) % len(p3)])
        return v

    nflip = 0
    for f in faces:
        p3 = np.array([bpos[r] for r in f])
        if np.dot(fnormal(p3), n) < 0:
            nflip += 1          # 面法向朝下 ⇒ 相对平面是"背面"
    # 真折叠：统一朝上后，投影有向面积仍为负
    nfold = 0
    for f in faces:
        p3 = np.array([bpos[r] for r in f])
        q = np.array([uv[r] for r in f])
        s = sarea(q)
        up = np.dot(fnormal(p3), n) >= 0
        if (s < 0) == up:
            nfold += 1

    dev = np.abs((P - c) @ n)
    diam = 2 * np.sqrt(np.max(((P - c) ** 2).sum(axis=1)))
    shallow = dev.max() / diam

    # 连通分量（按面共享节点）
    n2f = defaultdict(set)
    for i, f in enumerate(faces):
        for r in f:
            n2f[r].add(i)
    seen, comps = set(), []
    for i in range(len(faces)):
        if i in seen:
            continue
        st, comp = [i], set()
        while st:
            x = st.pop()
            if x in seen:
                continue
            seen.add(x)
            comp.add(x)
            for r in faces[x]:
                st.extend(n2f[r] - seen)
        comps.append(len(comp))
    comps.sort(reverse=True)

    ec = Counter()
    for f in faces:
        for i in range(len(f)):
            a, b = f[i], f[(i + 1) % len(f)]
            ec[(min(a, b), max(a, b))] += 1
    bnd = [k for k, v in ec.items() if v == 1]
    adj = defaultdict(set)
    for a, b in bnd:
        adj[a].add(b)
        adj[b].add(a)
    seen2, loops = set(), 0
    for a in list(adj):
        if a in seen2:
            continue
        loops += 1
        st = [a]
        while st:
            x = st.pop()
            if x in seen2:
                continue
            seen2.add(x)
            st.extend(adj[x] - seen2)

    print(f"{j:22s} {len(nodes):5d} {len(faces):5d} {len(comps):6d} "
          f"{str(comps[:4]):>18s} {shallow:7.3f} {nfold:6d} {loops:6d} "
          f"{100*len(bnd)/len(ec):7.1f}%")
print("=" * 96)
print("判据：连通块=1、真折叠=0、边界环=1、边界占比低 ⇒ 可参数化")
