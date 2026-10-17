"""A5: 从 .xplt(h5) 提取应力，与吴恺 3.74 MPa 量级对照。"""
from pathlib import Path
import numpy as np
import h5py

H5 = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio\ankle_tibiotalar.h5")

# pyfebio 的域标签有错位，按单元数映射回真实域名
BY_COUNT = {2415: "tibia_cort(hex8)", 6: "tibia_cort(tet4/楔拆分)",
            6651: "tibia_end_spon(松质)", 14482: "tibia_marrow(骨髓)",
            394: "talus_cort(hex8)", 1766: "talus_spon(松质)",
            458: "tibiotalar_tibia(软骨)", 613: "tibiotalar_talus(软骨)"}


def von_mises(s):
    sxx, syy, szz, sxy, syz, sxz = (s[:, i] for i in range(6))
    return np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2)
                   + 3.0 * (sxy ** 2 + syz ** 2 + sxz ** 2))


with h5py.File(H5, "r") as f:
    st = "13"          # 最后一个 state = 满载荷
    grp = f[f"states/{st}/element_data/stress"]
    print(f"=== state {st} 各域 von Mises 应力 (MPa) ===")
    print(f"{'域':26s} {'单元':>7s} {'min':>9s} {'中位':>9s} {'p99':>9s} {'max':>9s}")
    allvm = []
    for k in grp:
        arr = np.asarray(grp[k])
        vm = von_mises(arr)
        nm = BY_COUNT.get(arr.shape[0], k)
        allvm.append((nm, vm))
        print(f"{nm:26s} {arr.shape[0]:>7,} {vm.min():>9.4f} {np.median(vm):>9.4f} "
              f"{np.percentile(vm, 99):>9.4f} {vm.max():>9.4f}")

    flat = np.concatenate([v for _, v in allvm])
    print(f"\n全场 von Mises: 中位 {np.median(flat):.4f}  p99 {np.percentile(flat,99):.4f} "
          f"max {flat.max():.4f} MPa")

    # 最大 von Mises 出现在哪个域
    top = max(allvm, key=lambda kv: kv[1].max())
    print(f"最大应力出现在: {top[0]}")

    # 软骨应力（关节接触力学最关心）
    cart = [(n, v) for n, v in allvm if "软骨" in n]
    if cart:
        print(f"\n软骨最大 von Mises: " +
              " / ".join(f"{n} {v.max():.4f}" for n, v in cart))

    # 位移
    disp = np.asarray(f[f"states/{st}/node_data/displacement/all"])
    mag = np.linalg.norm(disp, axis=1)
    print(f"\n位移: max {mag.max():.4f} mm  中位 {np.median(mag):.4f} mm")
    print(f"  轴向(加载方向)分量 max {disp[:, 0].max():.4f} / min {disp[:, 0].min():.4f}")

print(f"""
=== 与吴恺 2012 对照 ===
  吴恺: 600 N 轴向载荷，软骨/骨接触区 von Mises ≈ 3.74 MPa
  本次: 顶面指定位移 {0.05} mm 压缩，量级见上
""")
