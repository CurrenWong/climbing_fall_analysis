"""查四边底行的**偏移节点配对顺序**：底 (b0,b1,b2,b3) 与 offset 哪种配对才得到正 det。

候选配对（对六面体求符号体积 det[b0,b1,b2,b3 | o_a,o_b,o_c,o_d]）：
  ① 同序      o0,o1,o2,o3
  ② 反序      o3,o2,o1,o0
  ③ 循环右移  o1,o2,o3,o0
  ④ 循环左移  o3,o0,o1,o2
"""
import numpy as np
from pathlib import Path

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
cg = np.load(W / "cartilage_mesh.npz", allow_pickle=True)
nid = cg["bone_node_ids"]
xyz = np.asarray(cg["bone_xyz"], float)
idx = {int(v): i for i, v in enumerate(nid)}


def hex_signed_vol(B, T):
    """六面体符号体积：拆成 6 个 tet（绕 0-6 对角线）。"""
    pts = np.vstack([B, T])          # 0..7
    HEX_T = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
             (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]
    s = 0.0
    for a, b, c, d in HEX_T:
        s += np.linalg.det(np.stack([pts[b] - pts[a], pts[c] - pts[a],
                                     pts[d] - pts[a]]))
    return s


for j in ("tibiotalar_tibia", "tibiotalar_talus"):
    E = cg[f"{j}_elems"]
    O = np.asarray(cg[f"{j}_offsets"], float)
    res = {k: [] for k in ("same", "rev", "rot1", "rot3")}
    for row in E:
        vals = [int(x) for x in row if int(x) != 0]
        if not vals or vals[0] <= 0:
            continue
        m = len(vals) // 2
        if m != 4:
            continue
        base = vals[:m]
        if any(b not in idx for b in base):
            continue
        B = np.array([xyz[idx[b]] for b in base])
        T = O[[-x - 1 for x in vals[m:]]]
        res["same"].append(hex_signed_vol(B, T[np.array([0, 1, 2, 3])]))
        res["rev"].append(hex_signed_vol(B, T[np.array([3, 2, 1, 0])]))
        res["rot1"].append(hex_signed_vol(B, T[np.array([1, 2, 3, 0])]))
        res["rot3"].append(hex_signed_vol(B, T[np.array([3, 0, 1, 2])]))
    print(f"\n{'='*74}\n{j}: 四边底行 {len(res['same'])} 个")
    print(f"  {'配对':>8} {'>0 个数':>9} {'<0 个数':>9} {'中位 |V|':>12}")
    for k in ("same", "rev", "rot1", "rot3"):
        a = np.array(res[k])
        print(f"  {k:>8} {(a>0).sum():9d} {(a<0).sum():9d} "
              f"{np.median(np.abs(a)):12.3f}")
    print("  ⇒ 只有「<0 个数 = 0」的那个配对是正确朝向")
