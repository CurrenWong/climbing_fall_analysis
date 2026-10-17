"""ATFL 取证：A/B 节点到底在解剖什么位置？为何相距 34 mm。"""
import numpy as np
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
OUT = ROOT / "temp/thums"
z = np.load(OUT / "thums_lowerlimb.npz", allow_pickle=True)
nid, xyz = z["node_ids"], z["node_xyz"]
i_of = {int(v): i for i, v in enumerate(nid)}

BONE_OF_PID = {
    81000700: "tibia_cort", 81000600: "tibia_end_spon", 81000601: "tibia_marrow",
    81000900: "fibula_cort", 81000800: "fibula_end_spon",
    81001100: "talus_cort", 81001000: "talus_spon",
    81001300: "calcaneus_cort", 81001200: "calcaneus_spon",
    81001400: "navicular_spon", 81001401: "navicular_cort",
}
BONE_NODES = {}
for pid, dom in BONE_OF_PID.items():
    s = set()
    for kind in ("solid", "shell"):
        k = f"{pid:08d}_{kind}"
        if k in z:
            s |= {int(v) for v in z[k].ravel()}
    BONE_NODES.setdefault(dom.split("_")[0], set()).update(s)
BONE_XYZ = {b: xyz[[i_of[n] for n in sorted(s)]] for b, s in BONE_NODES.items()}


def nearest_bone(P):
    """每个点最近的骨 + 距离。"""
    out = []
    for q in P:
        best = (None, 1e9)
        for b, Q in BONE_XYZ.items():
            d = float(np.sqrt(((Q - q) ** 2).sum(1).min()))
            if d < best[1]:
                best = (b, d)
        out.append(best)
    return out


for pid, nm, a_n, b_n in ((81100901, "ATFL", "talus", "fibula"),
                          (81100801, "CFL", "calcaneus", "fibula")):
    e = z[f"{pid:08d}_shell"]
    nodes = sorted({int(v) for v in e.ravel()})
    P = xyz[[i_of[n] for n in nodes]]
    print("=" * 100)
    print(f"=== {nm} (pid {pid})  壳单元 {len(e)} 个 / 节点 {len(nodes)} 个")
    # 精确延伸（min/max 投影）
    c = P - P.mean(0)
    _, _, vt = np.linalg.svd(c, full_matrices=False)
    for k in range(3):
        pr = c @ vt[k]
        print(f"  PCA 轴{k}: 延伸 {pr.max()-pr.min():>7.2f} mm")
    A = [n for n in nodes if n in BONE_NODES[a_n]]
    B = [n for n in nodes if n in BONE_NODES[b_n]]
    both = set(A) & set(B)
    A = [n for n in A if n not in both]
    B = [n for n in B if n not in both]
    rest = [n for n in nodes if n not in A and n not in B]
    print(f"  A({a_n})={len(A)}  B({b_n})={len(B)}  重合={len(both)}  未分配={len(rest)}")
    for lbl, grp in ((f"A={a_n}", A), (f"B={b_n}", B), ("未分配", rest)):
        if not grp:
            continue
        Q = xyz[[i_of[n] for n in grp]]
        print(f"    {lbl:<12} 质心={np.round(Q.mean(0),1)}  "
              f"包围盒对角={np.round(Q.max(0)-Q.min(0),1)}")
    print(f"  逐节点（按到 {b_n} 的距离升序）：")
    nb = nearest_bone(P)
    order = np.argsort([d for _, d in nb])
    print(f"    {'THUMS id':>10} {'坐标':<28} {'最近骨':<11} {'距':>6}  "
          f"{'∈A':>3} {'∈B':>3}")
    for k in order:
        n = nodes[k]
        b, d = nb[k]
        print(f"    {n:>10} {str(np.round(P[k],1)):<28} {b:<11} {d:>6.2f}  "
              f"{'✓' if n in A else '':>3} {'✓' if n in B else '':>3}")
