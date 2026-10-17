"""精确量化 `[ERR-20261006-003]` 修复带来的刚度变化（旧 PCA 轴 vs 新骨-骨轴）。"""
import io
import numpy as np
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
z = np.load(ROOT / "temp/thums/thums_lowerlimb.npz", allow_pickle=True)
nid, xyz = z["node_ids"], z["node_xyz"]
i_of = {int(v): i for i, v in enumerate(nid)}

BONE_OF_PID = {
    81000700: "tibia", 81000600: "tibia", 81000601: "tibia",
    81000900: "fibula", 81000800: "fibula",
    81001100: "talus", 81001000: "talus",
    81001300: "calcaneus", 81001200: "calcaneus",
    81001400: "navicular", 81001401: "navicular",
}
BONE_NODES = {}
for pid, dom in BONE_OF_PID.items():
    s = set()
    for kind in ("solid", "shell"):
        k = f"{pid:08d}_{kind}"
        if k in z:
            s |= {int(v) for v in z[k].ravel()}
    BONE_NODES.setdefault(dom, set()).update(s)

LIG_DEF = [
    (81100701, "PTFL", "talus", "fibula", 180.0, 2.0, "20~25"),
    (81100801, "CFL", "calcaneus", "fibula", 510.0, 2.0, "28.5"),
    (81100901, "ATFL", "talus", "fibula", 180.0, 2.0, "18~25"),
    (81101001, "PTCL", "talus", "calcaneus", 75.4, 2.0, "—"),
    (81101101, "SYND", "tibia", "fibula", 75.4, 2.0, "—"),
    (81101201, "MD1", "tibia", "talus", 180.0, 2.0, "—"),
    (81101202, "MD2", "tibia", "calcaneus", 510.0, 2.0, "—"),
    (81101203, "MD3", "tibia", "navicular", 510.0, 2.0, "—"),
    (81101204, "MD4", "tibia", "talus", 180.0, 2.0, "—"),
    (81101409, "SPRING", "calcaneus", "navicular", 75.4, 2.0, "—"),
]

print(f"{'韧带':<7} {'L_旧(PCA)':>10} {'L_新(骨-骨)':>11} {'L_真实杆长':>10} "
      f"{'A_旧':>8} {'A_新':>8} {'刚度变化':>9}  {'解剖':<8}")
print("─" * 100)
for pid, nm, a_n, b_n, ea, th, anat in LIG_DEF:
    e = z[f"{pid:08d}_shell"]
    nodes = sorted({int(v) for v in e.ravel()})
    P = xyz[[i_of[n] for n in nodes]]
    A = [n for n in nodes if n in BONE_NODES[a_n]]
    B = [n for n in nodes if n in BONE_NODES[b_n]]
    both = set(A) & set(B)
    A = [n for n in A if n not in both]
    B = [n for n in B if n not in both]
    c = P - P.mean(0)
    # 面积
    area = 0.0
    for q in e:
        p = xyz[[i_of[int(n)] for n in q]]
        area += 0.5 * np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0]))
        area += 0.5 * np.linalg.norm(np.cross(p[2] - p[0], p[3] - p[0]))
    # 旧：PCA 轴
    _, _, vt = np.linalg.svd(c, full_matrices=False)
    ax_old = vt[0]
    pr = c @ ax_old
    L_old = float(pr.max() - pr.min())
    # 新：骨-骨轴
    Pa = xyz[[i_of[n] for n in A]] if A else np.zeros((0, 3))
    Pb = xyz[[i_of[n] for n in B]] if B else np.zeros((0, 3))
    if len(Pa) and len(Pb):
        ax_new = (Pb.mean(0) - Pa.mean(0))
        ax_new /= np.linalg.norm(ax_new)
        pr2 = c @ ax_new
        L_new = float(pr2.max() - pr2.min())
        ang = float(np.degrees(np.arccos(np.clip(abs(ax_old @ ax_new), -1, 1))))
    else:
        L_new, ang = float("nan"), float("nan")
    A_old = area / L_old * th if L_old > 0 else 0.0
    A_new = area / L_new * th if L_new and L_new > 0 else 0.0
    ratio = A_new / A_old if A_old > 0 else float("nan")
    print(f"{nm:<7} {L_old:>10.2f} {L_new:>11.2f} {'':>10} "
          f"{A_old:>8.2f} {A_new:>8.2f} {ratio:>9.2f}× {anat:<8}"
          f"   轴夹角 {ang:.1f}°")

print("\n★ 刚度 ∝ A_cross/L 的净变化要看 EA：FEBio 杆用 (EA/L₀) 作轴向刚度，")
print("  其中 L₀ = 杆的实际几何长度（配对决定，未变），A 由上面的 A_each 决定。")
print("  ⇒ 轴修正在 **PTCL/SYND**（夹角 84~89°）上效果最大。")
