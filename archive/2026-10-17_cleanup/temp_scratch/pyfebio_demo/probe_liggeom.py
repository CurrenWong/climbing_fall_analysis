"""韧带几何诊断：找出 ATFL 为何 36.99 mm（解剖 15~20）。

嫌疑：`ligaments.py` 用**整个点云的主成分最大奇异向量**当"纤维轴"：
    c = P - P.mean(0); _, _, vt = np.linalg.svd(c); ax = vt[0]
      ⇒ 对**宽而短**的韧带（ATFL 正是"宽带"型），最大延伸可能**横跨宽度**而非沿纤维
      ⇒ ① 轴错了 ② 沿轴配对会把 A/B 节点**横向错配** ⇒ 杆长≠真纤维长
本脚本把证据摆出来：三条 PCA 延伸、A/B 质心距、主成分轴 vs (B−A) 方向夹角、杆长分布。
"""
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
LIG_DEF = [
    (81100701, "PTFL", "talus", "fibula", "20~25"),
    (81100801, "CFL", "calcaneus", "fibula", "25~32 (28.5)"),
    (81100901, "ATFL", "talus", "fibula", "**15~20**"),
    (81101001, "PTCL", "talus", "calcaneus", "—"),
    (81101101, "SYND", "tibia", "fibula", "—"),
    (81101201, "MD1", "tibia", "talus", "—"),
    (81101202, "MD2", "tibia", "calcaneus", "—"),
    (81101203, "MD3", "tibia", "navicular", "—"),
    (81101204, "MD4", "tibia", "talus", "—"),
    (81101409, "SPRING", "calcaneus", "navicular", "—"),
]
BONE_NODES = {}
for pid, dom in BONE_OF_PID.items():
    s = set()
    for kind in ("solid", "shell"):
        k = f"{pid:08d}_{kind}"
        if k in z:
            s |= {int(v) for v in z[k].ravel()}
    BONE_NODES.setdefault(dom.split("_")[0], set()).update(s)

print(f"{'韧带':<7} {'节点':>4} {'PCA延伸(mm) 大/中/小':<26} "
      f"{'|A-B|':>7} {'夹角°':>6} {'杆长 min/mean/max':<22} {'解剖':<12} 判定")
print("─" * 132)
for pid, nm, a_name, b_name, anat in LIG_DEF:
    e = z[f"{pid:08d}_shell"]
    nodes = sorted({int(v) for v in e.ravel()})
    P = xyz[[i_of[n] for n in nodes]]
    A = [n for n in nodes if n in BONE_NODES[a_name]]
    B = [n for n in nodes if n in BONE_NODES[b_name]]
    both = set(A) & set(B)
    A = [n for n in A if n not in both]
    B = [n for n in B if n not in both]

    c = P - P.mean(0)
    sv = np.linalg.svd(c, compute_uv=False)          # 3 个奇异值 ≈ 延伸量
    ext = 2 * sv / np.sqrt(len(P))                   # 粗略"半宽"×2
    Pa = xyz[[i_of[n] for n in A]] if A else np.zeros((0, 3))
    Pb = xyz[[i_of[n] for n in B]] if B else np.zeros((0, 3))
    ca = Pa.mean(0) if len(Pa) else np.zeros(3)
    cb = Pb.mean(0) if len(Pb) else np.zeros(3)
    dAB = float(np.linalg.norm(cb - ca)) if len(Pa) and len(Pb) else float("nan")
    _, _, vt = np.linalg.svd(c, full_matrices=False)
    ax_pca = vt[0] / np.linalg.norm(vt[0])
    ax_bone = (cb - ca)
    ax_bone = ax_bone / max(np.linalg.norm(ax_bone), 1e-12)
    ang = float(np.degrees(np.arccos(np.clip(abs(ax_pca @ ax_bone), -1, 1))))
    # 现有配对法的杆长
    pr = c @ ax_pca
    if len(Pa) and len(Pb):
        tA = np.array([(xyz[i_of[n]] - P.mean(0)) @ ax_pca for n in A])
        tB = np.array([(xyz[i_of[n]] - P.mean(0)) @ ax_pca for n in B])
        lo, hi = pr.min(), pr.max()
        tA = (tA - lo) / (hi - lo + 1e-12)
        tB = (tB - lo) / (hi - lo + 1e-12)
        Ls = []
        for n, t in zip(A, tA):
            j = int(np.argmin(np.abs(tB - t)))
            Ls.append(float(np.linalg.norm(xyz[i_of[B[j]]] - xyz[i_of[n]])))
        Ls = np.array(Ls)
        lstr = f"{Ls.min():.1f}/{Ls.mean():.1f}/{Ls.max():.1f}"
    else:
        lstr = "—"
    # 判定：主成分轴是否与骨-骨方向差很多
    verdict = ""
    if dAB == dAB and ang > 25:
        verdict = f"⚠ 主成分轴与骨-骨方向差 {ang:.0f}° ⇒ 轴选错"
    elif nm == "ATFL" and Ls.mean() > 25:
        verdict = "⚠ 杆长偏长（配对错）"
    print(f"{nm:<7} {len(nodes):>4} {ext[0]:>7.1f}/{ext[1]:>6.1f}/{ext[2]:>6.1f}      "
          f"{dAB:>7.2f} {ang:>6.1f} {lstr:<22} {anat:<12} {verdict}")
