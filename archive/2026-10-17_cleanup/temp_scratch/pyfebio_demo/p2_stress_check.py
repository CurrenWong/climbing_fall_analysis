"""P2d：查 σ_max=104.85 是真奇点还是读数伪影（原始 Voigt 分量 + 跨载荷一致性）。"""
from pathlib import Path

import numpy as np
import h5py

R = Path(r"D:\Project\climbing_fall_analysis")
H5 = R / "temp/thums/febio_lig/ankle_lig_WK_stress.h5"


def vm(s):
    s = np.asarray(s, float)
    xx, yy, zz, xy, yz, xz = (s[..., i] for i in range(6))
    return np.sqrt(0.5*((xx-yy)**2 + (yy-zz)**2 + (zz-xx)**2)
                   + 3.0*(xy**2 + yz**2 + xz**2))


with h5py.File(H5, "r") as f:
    sids = sorted(f["states"].keys(), key=int)
    sg = f["states"][sids[-1]]["element_data"]["stress"]
    for d in sorted(sg.keys()):
        a = np.asarray(sg[d][:], float)
        if "tibiotalar" not in d:
            continue
        print("=" * 90)
        print(f"域 {d}   shape={a.shape}  dtype={a.dtype}")
        if a.ndim == 2 and a.shape[1] == 6:
            v = vm(a)
            print(f"  vm: min {v.min():.3f} p50 {np.median(v):.3f} "
                  f"p99 {np.percentile(v,99):.3f} max {v.max():.3f}")
            j = int(np.argmax(v))
            print(f"  最大单元 #{j}: Voigt=[{a[j,0]:.2f} {a[j,1]:.2f} "
                  f"{a[j,2]:.2f} {a[j,3]:.2f} {a[j,4]:.2f} {a[j,5]:.2f}]  vm={v[j]:.2f}")
            top = np.argsort(v)[-8:][::-1]
            print("  前 8 大单元:")
            for k in top:
                print(f"    #{k:<5} vm={v[k]:>9.2f}  主对角="
                      f"[{a[k,0]:>8.2f} {a[k,1]:>8.2f} {a[k,2]:>8.2f}]")
            print(f"  各分量 max: xx {np.abs(a[:,0]).max():.2f} "
                  f"yy {np.abs(a[:,1]).max():.2f} zz {np.abs(a[:,2]).max():.2f} "
                  f"xy {np.abs(a[:,3]).max():.2f} yz {np.abs(a[:,4]).max():.2f} "
                  f"xz {np.abs(a[:,5]).max():.2f}")
        else:
            print(f"  非 (n,6) 形状 ⇒ 口径不是 Voigt！")

print("\n" + "=" * 90)
print("对照：其他域的同口径")
print("=" * 90)
with h5py.File(H5, "r") as f:
    sids = sorted(f["states"].keys(), key=int)
    sg = f["states"][sids[-1]]["element_data"]["stress"]
    for d in sorted(sg.keys()):
        a = np.asarray(sg[d][:], float)
        if a.ndim != 2 or a.shape[1] != 6:
            continue
        v = vm(a)
        print(f"  {d:<34} n={len(v):>5}  max {v.max():>9.3f}  "
              f"p99 {np.percentile(v,99):>8.3f}  p50 {np.median(v):>8.3f}")

print("\n" + "=" * 90)
print("跨 state 一致性（同一域的 max 随载荷单调吗）")
print("=" * 90)
with h5py.File(H5, "r") as f:
    sids = sorted(f["states"].keys(), key=int)
    for d in ("tibiotalar_talus__penta6", "tibiotalar_tibia__penta6"):
        seq = []
        for s in sids:
            g = f["states"][s]["element_data"]["stress"]
            if d in g:
                v = vm(np.asarray(g[d][:], float))
                seq.append((s, float(v.max()), float(np.percentile(v, 99))))
        print(f"  {d}")
        for s, mx, p99 in seq:
            print(f"    state {s:>3}: max {mx:>10.3f}  p99 {p99:>8.3f}")
