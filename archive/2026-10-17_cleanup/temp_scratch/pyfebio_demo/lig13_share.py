"""决定性检查：韧带节点是否与骨节点共用 ID？顺带算每条韧带的几何刚度先验。"""
import numpy as np
from pathlib import Path

T = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
z = np.load(T / "thums_lowerlimb.npz", allow_pickle=True)
nid, xyz = z["node_ids"], z["node_xyz"]
i_of = {int(v): i for i, v in enumerate(nid)}

BONE_KINDS = {}
for pid in (81000700, 81000600, 81000900, 81000800, 81001100, 81001000,
            81001300, 81001200, 81001400, 81001401):
    for kind in ("solid", "shell"):
        k = f"{pid:08d}_{kind}"
        if k in z:
            BONE_KINDS[k] = z[k]

bone_nodes = set()
for e in BONE_KINDS.values():
    bone_nodes |= {int(v) for v in e.ravel()}
print(f"骨节点合计 {len(bone_nodes):,}")

LIG = [(81100701, "PTFL 距腓后", "shell"), (81100801, "CFL 跟腓", "shell"),
       (81100901, "ATFL 距腓前", "shell"),
       (81101001, "PTCL 距跟后", "shell"),
       (81101101, "下胫腓", "shell"),
       (81101201, "内侧1 胫距后深", "shell"),
       (81101202, "内侧2 胫跟", "shell"),
       (81101203, "内侧3 胫舟", "shell"),
       (81101204, "内侧4 胫距前深", "shell"),
       (81101401, "跖长", "shell")]

EA = {81100701: 180, 81100801: 510, 81100901: 180, 81101001: 75.4,
      81101101: 75.4, 81101201: 180, 81101202: 510, 81101203: 510,
      81101204: 180, 81101401: 75.4}

print(f"\n{'pid':>9} {'名称':16s} {'单元':>4} {'节点':>4} {'共用骨节点':>10} "
      f"{'未共用':>7} {'跨度mm':>7} {'面积mm²':>8} {'EA':>6}")
print("-" * 92)
rows = []
for pid, nm, kind in LIG:
    k = f"{pid:08d}_{kind}"
    e = z[k]
    nodes = {int(n) for n in e.ravel()}
    shared = nodes & bone_nodes
    # 几何：用质心分解面积 + 最大端距作为长度
    P = xyz[[i_of[n] for n in sorted(nodes)]]
    # 面积：把每个 quad 拆 2 三角
    area = 0.0
    for q in e:
        p = xyz[[i_of[int(n)] for n in q]]
        area += 0.5 * np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0]))
        area += 0.5 * np.linalg.norm(np.cross(p[2] - p[0], p[3] - p[0]))
    # 长度：PCA 主轴范围
    c = P - P.mean(0)
    u, s, vt = np.linalg.svd(c, full_matrices=False)
    proj = c @ vt[0]
    length = proj.max() - proj.min()
    rows.append((pid, nm, len(e), len(nodes), len(shared), len(nodes) - len(shared),
                 length, area, EA[pid]))
    print(f"{pid:>9} {nm:16s} {len(e):>4} {len(nodes):>4} {len(shared):>10} "
          f"{len(nodes)-len(shared):>7} {length:>7.1f} {area:>8.2f} {EA[pid]:>6}")

print("\n\n=== 几何刚度先验 k = EA*A/L（N/mm），用于预判排序 ===")
print("（仅数量级参考；真实排序还要看加载方向投影）\n")
sc = []
for pid, nm, ne, nn, sh, un, L, A, ea in rows:
    k = ea * A / L if L > 0 else 0
    sc.append((k, pid, nm, ea, A, L))
for k, pid, nm, ea, A, L in sorted(sc, reverse=True):
    print(f"  {k:9.1f} N/mm   {nm:16s} EA={ea:5} A={A:6.2f} L={L:5.1f}")

print("\n=== 吴恺要的四条 ===")
want = {81101202: "胫跟", 81101203: "胫舟", 81100801: "跟腓",
        81101101: "下胫腓前"}
order = sorted([r for r in sc if r[1] in want], reverse=True)
print(f"  吴恺排序: 胫跟 > 胫舟 > 跟腓 > 下胫腓前")
print(f"  几何先验: " + " > ".join(f"{want[p]}" for k, p, nm, *rest in order))
