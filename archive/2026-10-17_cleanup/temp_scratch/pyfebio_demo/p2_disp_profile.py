"""P2e（最关键）：1.6mm 指定位移的位移剖面 —— 柔度到底在哪。"""
from pathlib import Path

import numpy as np
import h5py

R = Path(r"D:\Project\climbing_fall_analysis")
H5 = R / "temp/thums/febio_lig/ankle_lig_WK_stress.h5"

with h5py.File(H5, "r") as f:
    nd = f["meshes/0/nodes"][:]
    X = np.column_stack([nd["x"], nd["y"], nd["z"]]).astype(float)
    sids = sorted(f["states"].keys(), key=int)
    D = np.asarray(f["states"][sids[-1]]["node_data"]["displacement"]["all"][:],
                   float)
    # 域元素集
    doms = {}
    for dn in f["meshes/0/domains"].keys():
        try:
            es = f["meshes/0/domains"][dn]["ElementSet"]
            doms[dn] = np.asarray(es[:], int).ravel()
        except Exception:
            pass

# 末态坐标
Xd = X + D
dz = D[:, 2]          # 轴向（ALIGN_AXIS 后 z = 加载轴）

print("=" * 100)
print("位移剖面：按参考 z 分箱，看轴向位移 dz 随高度怎么长")
print("=" * 100)
z = X[:, 2]
lo, hi = np.percentile(z, 0.5), np.percentile(z, 99.5)
print(f"z 范围 {z.min():.2f} .. {z.max():.2f} mm")
edges = np.linspace(z.min(), z.max(), 22)
print(f"{'z 区间':>22}{'节点':>7}{'dz 均值':>11}{'dz 中位':>11}{'|u|均值':>11}")
for i in range(len(edges) - 1):
    m = (z >= edges[i]) & (z < edges[i + 1])
    if m.sum() < 5:
        continue
    print(f"[{edges[i]:>8.2f},{edges[i+1]:>8.2f}){m.sum():>7}"
          f"{dz[m].mean():>11.4f}{np.median(dz[m]):>11.4f}"
          f"{np.linalg.norm(D[m], axis=1).mean():>11.4f}")

print("\n" + "=" * 100)
print("关键三处的轴向位移")
print("=" * 100)


def dom_of(name):
    for k, v in doms.items():
        if k.startswith(name):
            return v
    return None


for nm in ("tibia_cort__hex8", "talus_cort__hex8", "tibiotalar_tibia__penta6",
           "tibiotalar_talus__penta6"):
    e = dom_of(nm)
    if e is None:
        continue
    with h5py.File(H5, "r") as f:
        conn = np.asarray(f["meshes/0/domains"][nm]["ElementSet"][:], int).ravel()
    # 用元素集对应的节点：从 elementsets 取
    try:
        with h5py.File(H5, "r") as f:
            es = f["meshes/0/elementsets"]
            key = [k for k in es.keys() if nm.split("__")[0] in k]
    except Exception:
        pass

# 直接用全体节点的 z 分布 + 位移（最简）
print(f"\n  胫骨顶面（z 最大 2%）: dz = {dz[z >= np.percentile(z,98)].mean():.4f} mm")
print(f"  距骨底（z 最小 2%）  : dz = {dz[z <= np.percentile(z,2)].mean():.4f} mm")

print("\n" + "=" * 100)
print("与理论刚度对照")
print("=" * 100)
F = 592.602
d_top = dz[z >= np.percentile(z, 98)].mean()
print(f"  实测总位移（顶面 dz）= {d_top:.4f} mm   F = {F:.1f} N")
print(f"  ⇒ k_total = {F/d_top:.0f} N/mm")
print()
print("  若柔度全在胫骨（E=12000 MPa, A≈300 mm², L=60 mm）：")
k_t = 12000 * 300 / 60
print(f"     k_tibia = {k_t:.0f} N/mm ⇒ {F/k_t:.4f} mm（远小于实测）")
print("  若柔度全在关节软骨（E=12 MPa, A≈300 mm², t=2×1.3 mm）：")
k_j = 12 * 300 / 2.6
print(f"     k_joint = {k_j:.0f} N/mm ⇒ {F/k_j:.4f} mm")
print(f"\n  ⇒ 实测 {d_top:.3f} mm 比这两者都大 ⇒ 柔度来自**别处**（韧带？）")
