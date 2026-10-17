"""P2f：按域定位柔度（每个域的 dz 跨度 / 轴向应变 / 隐含刚度）。"""
from pathlib import Path

import numpy as np
import h5py

R = Path(r"D:\Project\climbing_fall_analysis")
H5 = R / "temp/thums/febio_lig/ankle_lig_WK_stress.h5"
F = 592.602

with h5py.File(H5, "r") as f:
    nd = f["meshes/0/nodes"][:]
    X = np.column_stack([nd["x"], nd["y"], nd["z"]]).astype(float)
    sids = sorted(f["states"].keys(), key=int)
    D = np.asarray(f["states"][sids[-1]]["node_data"]["displacement"]["all"][:],
                   float)
    # 网格域的节点集
    nodesets = {}
    if "nodesets" in f["meshes/0"]:
        for k in f["meshes/0/nodesets"].keys():
            nodesets[k] = np.asarray(f["meshes/0/nodesets"][k][:], int).ravel()
    domains = list(f["meshes/0/domains"].keys())
    print("domains:", domains)
    print("nodesets:", list(nodesets.keys())[:20])

dz = D[:, 2]
zr = X[:, 2] - X[:, 2].min()

print("\n" + "=" * 100)
print("逐域：z 跨度 / dz 跨度 / 轴向应变估计")
print("=" * 100)
print(f"{'域':<32}{'节点':>7}{'z跨度':>10}{'dz跨度':>11}{'dz均值':>10}"
      f"{'|u|max':>10}")
for dn in domains:
    try:
        es = f["meshes/0/domains"][dn]["ElementSet"]
        el = np.asarray(es[:], int).ravel()
    except Exception:
        continue
    # 元素→节点：从 elementsets 找 connectivity
    try:
        with h5py.File(H5, "r") as ff:
            et = ff["meshes/0/elementsets"][dn]
            conn = np.asarray(et[:])
    except Exception:
        conn = None
    nodes = None
    if conn is not None and conn.ndim == 2:
        nodes = np.unique(conn.ravel()) - 1
    if nodes is None or len(nodes) < 3:
        continue
    zz, dd = X[nodes, 2], dz[nodes]
    print(f"{dn:<32}{len(nodes):>7}{zz.max()-zz.min():>10.2f}"
          f"{dd.max()-dd.min():>11.4f}{dd.mean():>10.4f}"
          f"{np.linalg.norm(D[nodes],axis=1).max():>10.4f}")

print("\n" + "=" * 100)
print("韧带域（__unnamed_*）：两端 dz 差 = 拉伸量")
print("=" * 100)
with h5py.File(H5, "r") as f:
    for dn in domains:
        if not dn.startswith("__unnamed"):
            continue
        try:
            conn = np.asarray(f["meshes/0/elementsets"][dn][:])
        except Exception:
            continue
        nodes = np.unique(conn.ravel()) - 1
        if len(nodes) < 2:
            continue
        zz, dd = X[nodes, 2], dz[nodes]
        # 两端 = z 最大/最小的节点
        hi = nodes[np.argmax(zz)]
        lo = nodes[np.argmin(zz)]
        ztop, zbot = X[hi, 2], X[lo, 2]
        dtop, dbot = dz[hi], dz[lo]
        L = ztop - zbot
        stretch = dtop - dbot
        strain = stretch / L if L else float("nan")
        print(f"{dn:<18} 节点{len(nodes):>4}  z跨度 {L:>7.2f} mm  "
              f"Δdz {stretch:>8.4f} mm  εz {strain*100:>7.2f}%  "
              f"|u|max {np.linalg.norm(D[nodes],axis=1).max():.3f}")

print("\n" + "=" * 100)
print("👉 判定：哪个域贡献了 1.6mm 中的绝大部分")
print("=" * 100)
print(f"  总指定位移 1.6000 mm，关节只闭 0.0142 mm ⇒ 其余 "
      f"{1.6-0.0142:.4f} mm 在别处")
