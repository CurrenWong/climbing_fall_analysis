"""探 cartilage_mesh.npz 的精确结构 + 底面面型构成 + 现有厚度场。"""
import numpy as np
from pathlib import Path
from collections import Counter

OUT = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
NPZ = OUT / "cartilage_mesh.npz"
cg = np.load(NPZ)
print(f"文件 {NPZ.name}  keys({len(cg.files)}):")
for k in cg.files:
    a = cg[k]
    print(f"   {k:34s} shape={str(a.shape):14s} dtype={a.dtype}")

db = np.load(OUT / "thums_lowerlimb.npz")
nid, nxyz = db["node_ids"], db["node_xyz"]
idx = {int(v): i for i, v in enumerate(nid)}

CARTS = [k.replace("_elems", "") for k in cg.files if k.endswith("_elems")]
print(f"\n关节层 {len(CARTS)}: {CARTS}\n")
print(f"{'关节':24s} {'格数':>6s} {'tri':>5s} {'quad':>5s} | "
      f"{'底节点':>7s} {'偏移节点':>8s} {'外表面':>7s} | "
      f"{'厚min':>7s} {'厚中位':>7s} {'厚max':>7s}  单位mm")
for j in CARTS:
    e = cg[f"{j}_elems"]
    nk = f"{j}_nodes" if f"{j}_nodes" in cg.files else f"{j}_offsets"
    nd = cg[nk]
    fo = cg[f"{j}_is_offset"] if f"{j}_is_offset" in cg.files else None
    ntri = nquad = 0
    th = []
    for row in e:
        vals = [int(x) for x in row if int(x) != 0]
        if not vals:
            continue
        m = len(vals) // 2
        if m == 3:
            ntri += 1
        elif m == 4:
            nquad += 1
        base, offs = vals[:m], vals[m:]
        for b, o in zip(base, offs):
            if b > 0 and o < 0:
                k = -o - 1
                if k < len(nd):
                    th.append(float(np.linalg.norm(nd[k] - nxyz[idx[b]])))
    th = np.array(th)
    nout = int(fo.sum()) if fo is not None else len(nd)
    print(f"{j:24s} {ntri+nquad:6d} {ntri:5d} {nquad:5d} | "
          f"{len(nd):7d} {len(nd):8d} {nout:7d} | "
          f"{th.min():7.3f} {np.median(th):7.3f} {th.max():7.3f}")

print("\n★ 判据：tri 占比高说明底面不是全 quad（需对 tri 也做 1→4）")
print("★ 厚度中位应≈1.5 mm（Anderson 2006）")
