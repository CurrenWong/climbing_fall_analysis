"""从 THUMS 壳韧带生成「杆单元束」参数（对标吴恺：只承受拉应力的单轴梁单元）。

输出 temp/thums/ligaments.npz：
  - 每条韧带的节点（THUMS id + 坐标）
  - 两端附着节点（按骨分组）
  - 杆束配对 (a_id, b_id)  以及每根杆的截面积
  - 材料 EA / 壳厚

★ 截面积口径（易错）：
   THUMS 壳是**单层膜**，我算出的 294 mm² 是**膜面积 = 宽×长**，
   不是垂直于纤维的截面。真实截面 = 宽 × 厚 = (膜面积/长度) × 厚度。
   ⇒ 单根杆截面 = 该项 / 杆数。
"""
import json
import numpy as np
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
OUT = ROOT / "temp/thums"
z = np.load(OUT / "thums_lowerlimb.npz", allow_pickle=True)
nid, xyz = z["node_ids"], z["node_xyz"]
i_of = {int(v): i for i, v in enumerate(nid)}

# ── 骨部件（pid -> 域键, 名称）──
BONE_OF_PID = {
    81000700: "tibia_cort", 81000600: "tibia_end_spon",
    81000601: "tibia_marrow",
    81000900: "fibula_cort", 81000800: "fibula_end_spon",
    81001100: "talus_cort", 81001000: "talus_spon",
    81001300: "calcaneus_cort", 81001200: "calcaneus_spon",
    81001400: "navicular_spon", 81001401: "navicular_cort",
}

# ── 8 组踝韧带（吴恺）+ 弹簧韧带（稳定舟骨用）──
LIG_DEF = [
    # pid, 名称, EA(MPa), 厚度(mm)
    (81100701, "PTFL_talofibular_posterior", 180.0, 2.0),
    (81100801, "CFL_calcaneofibular", 510.0, 2.0),
    (81100901, "ATFL_talofibular_anterior", 180.0, 2.0),
    (81101001, "PTCL_talocalcaneal_posterior", 75.4, 2.0),
    (81101101, "SYND_tibiofibular", 75.4, 2.0),
    (81101201, "MD1_tibiotalar_post", 180.0, 2.0),
    (81101202, "MD2_tibiocalcaneal", 510.0, 2.0),
    (81101203, "MD3_tibionavicular", 510.0, 2.0),
    (81101204, "MD4_tibiotalar_ant", 180.0, 2.0),
    (81101409, "SPRING_calcaneonavicular", 75.4, 2.0),
]

# 各韧带允许的两侧骨（顺序 = A侧, B侧）
BONES_OF = {
    81100701: ("talus", "fibula"),
    81100801: ("calcaneus", "fibula"),
    81100901: ("talus", "fibula"),
    81101001: ("talus", "calcaneus"),
    81101101: ("tibia", "fibula"),
    81101201: ("tibia", "talus"),
    81101202: ("tibia", "calcaneus"),
    81101203: ("tibia", "navicular"),
    81101204: ("tibia", "talus"),
    81101409: ("calcaneus", "navicular"),
}

# 骨 -> 该骨的 THUMS 节点集合（仅用 solid+shell 的单元节点）
BONE_NODES = {}
for pid, dom in BONE_OF_PID.items():
    s = set()
    for kind in ("solid", "shell"):
        k = f"{pid:08d}_{kind}"
        if k in z:
            s |= {int(v) for v in z[k].ravel()}
    BONE_NODES.setdefault(dom.split("_")[0], set()).update(s)

print("=== 骨节点集合 ===")
for b, s in sorted(BONE_NODES.items()):
    print(f"  {b:12s} {len(s):>6,}")

records = []
print(f"\n{'pid':>9} {'名称':28s} {'节点':>4} {'A侧':>4} {'B侧':>4} "
      f"{'杆数':>4} {'膜面积':>8} {'长度':>6} {'截面':>7} {'每杆A':>7}")
print("-" * 100)
for pid, nm, ea, th in LIG_DEF:
    key = f"{pid:08d}_shell"
    e = z[key]
    nodes = sorted({int(v) for v in e.ravel()})
    P = xyz[[i_of[n] for n in nodes]]

    a_name, b_name = BONES_OF[pid]
    A = [n for n in nodes if n in BONE_NODES[a_name]]
    B = [n for n in nodes if n in BONE_NODES[b_name]]
    # 去掉同属两侧的
    both = set(A) & set(B)
    A = [n for n in A if n not in both]
    B = [n for n in B if n not in both]

    # 膜面积
    area = 0.0
    for q in e:
        p = xyz[[i_of[int(n)] for n in q]]
        area += 0.5 * np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0]))
        area += 0.5 * np.linalg.norm(np.cross(p[2] - p[0], p[3] - p[0]))

    # ── 纤维方向 = **骨-骨方向**（A 簇质心 → B 簇质心）★ 2026-10-06 修
    #   `[ERR-20261006-003]`：原实现用**整片点云的主成分最大奇异向量**当纤维轴。
    #   实测该轴对**宽而短**的韧带完全错（与骨-骨方向几乎垂直）：**PTCL 89.3°、SYND 84.0°**，
    #   于是 L 从真值 18.0/32.4 被抬到 50.7/44.5 ⇒ `A_cross = 面积/L×厚` 偏小
    #   ⇒ **轴向刚度偏低 PTCL 2.81× / SYND 1.37×**（已量化）。
    #   ⚠️ 反过来，长条型韧带本就不受影响：**CFL 轴夹角仅 0.7°、ATFL 2.1°、
    #   MD1~MD4 1.3~4.8° ⇒ 刚度变化 0.98~1.09×**（不要误以为"CFL 也被低估"）。
    #   韧带的纤维方向**按定义**是两端附着点连线，故用骨-骨方向，不用点云形状。
    Pa = xyz[[i_of[n] for n in A]] if A else np.zeros((0, 3))
    Pb = xyz[[i_of[n] for n in B]] if B else np.zeros((0, 3))
    ca = Pa.mean(0) if len(Pa) else np.zeros(3)
    cb = Pb.mean(0) if len(Pb) else np.zeros(3)
    c = P - P.mean(0)
    d_ab = float(np.linalg.norm(cb - ca)) if len(Pa) and len(Pb) else 0.0
    if d_ab > 1e-9:
        ax = (cb - ca) / d_ab                 # ★ 骨-骨方向（A→B）
        ax_src = "骨-骨"
    else:                                     # 回退：没有两端信息时只能用形状
        _, _, vt = np.linalg.svd(c, full_matrices=False)
        ax = vt[0] / np.linalg.norm(vt[0])
        ax_src = "PCA(回退)"
    pr = c @ ax
    L = float(pr.max() - pr.min())

    # 截面 = 宽 × 厚 = (膜面积/长) × 厚
    A_cross = area / L * th if L > 0 else 0.0

    # ── 杆束配对：**1:1 贪心最近邻**（3D 距离）★ 2026-10-06 修
    #   `[ERR-20261006-004]`：原实现「按沿轴归一化位置最近」在 ATFL 上退化成
    #   **4 根杆全部汇聚到同一个腓骨节点**（实测 `lig_pairs` = [[68,1497],[444,1497],
    #   [445,1497],[505,1497]]）⇒ **足迹被压成一个点** ⇒ 无法抵抗绕该点的转动，
    #   韧带对关节的转动约束被低估。
    #   改为：按 3D 距离升序贪心取用未被占用的端点 ⇒ 1:1 分布、足迹被分辨。
    trusses = []
    if len(Pa) and len(Pb):
        nA, nB = len(Pa), len(Pb)
        cand = sorted((float(np.linalg.norm(Pb[j] - Pa[i])), i, j)
                      for i in range(nA) for j in range(nB))
        usedA, usedB = set(), set()
        for d, i, j in cand:                   # ① 先 1:1 ⇒ 两侧足迹都被覆盖
            if i in usedA or j in usedB:
                continue
            usedA.add(i)
            usedB.add(j)
            trusses.append((A[i], B[j]))
        for i in range(nA):                    # ② 剩余 A 各接其最近 B（可复用）
            if i in usedA:
                continue
            j = int(np.argmin([np.linalg.norm(Pb[k] - Pa[i])
                               for k in range(nB)]))
            trusses.append((A[i], B[j]))
        for j in range(nB):                    # ③ 对称：剩余 B 各接其最近 A
            if j in usedB:
                continue
            i = int(np.argmin([np.linalg.norm(Pb[j] - Pa[k])
                               for k in range(nA)]))
            trusses.append((A[i], B[j]))
    ntr = len(trusses)
    A_each = A_cross / ntr if ntr else 0.0
    L_tru = np.array([float(np.linalg.norm(xyz[i_of[b]] - xyz[i_of[a]]))
                      for a, b in trusses]) if trusses else np.zeros(0)
    # ★ 自检（信息性，不是告警）：
    #   `L` = 膜沿**纤维轴**的延伸 = 面积公式 `A_cross = 面积/L×厚` 的正确分母（膜可能弯曲，
    #   所以直杆长会短于 L，属正常）。`[ERR-20261006-003]` 的真实影响只落在
    #   **主成分轴与骨-骨方向差很多**的韧带上（PTCL 89°、SYND 84°）。
    lmed = float(np.median(L_tru)) if len(L_tru) else float("nan")
    nA_cov = len({a for a, _ in trusses})
    nB_cov = len({b for _, b in trusses})

    print(f"{pid:>9} {nm:28s} {len(nodes):>4} {len(A):>4} {len(B):>4} "
          f"{ntr:>4} {area:>8.1f} {L:>6.1f} {A_cross:>7.2f} {A_each:>7.3f} "
          f"{ax_src:>7} {lmed:>7.1f}  覆盖A {nA_cov}/{len(A)} B {nB_cov}/{len(B)}")

    records.append(dict(pid=pid, name=nm, ea=ea, thick=th,
                        nodes=nodes, pid_key=key,
                        a_bone=a_name, b_bone=b_name,
                        A=A, B=B, trusses=trusses,
                        area=area, L=L, A_cross=A_cross, A_each=A_each))

# ── 打包 ──
save = {}
save["n_lig"] = np.array(len(records))
save["pids"] = np.array([r["pid"] for r in records])
save["names"] = np.array([r["name"] for r in records])
save["ea"] = np.array([r["ea"] for r in records])
save["thick"] = np.array([r["thick"] for r in records])
save["area_membrane"] = np.array([r["area"] for r in records])
save["length"] = np.array([r["L"] for r in records])
save["area_cross"] = np.array([r["A_cross"] for r in records])
save["a_bone"] = np.array([r["a_bone"] for r in records])
save["b_bone"] = np.array([r["b_bone"] for r in records])

mx = max(len(r["trusses"]) for r in records)
TR = -np.ones((len(records), mx, 2), dtype=np.int64)
NP_ = -np.ones((len(records), mx), dtype=np.float64)
NTR = np.zeros(len(records), dtype=np.int64)
for i, r in enumerate(records):
    NTR[i] = len(r["trusses"])
    NP_[i] = r["A_each"]
    for j, (a, b) in enumerate(r["trusses"]):
        TR[i, j] = (a, b)
save["truss"] = TR
save["truss_area"] = NP_
save["n_truss"] = NTR

# 所有韧带节点坐标
allv = sorted({n for r in records for n in r["nodes"]})
v_of = {n: i for i, n in enumerate(allv)}
save["v_thums_id"] = np.array(allv, dtype=np.int64)
save["v_xyz"] = xyz[[i_of[n] for n in allv]]
mxv = max(len(r["nodes"]) for r in records)
VI = -np.ones((len(records), mxv), dtype=np.int64)
NV = np.zeros(len(records), dtype=np.int64)
for i, r in enumerate(records):
    NV[i] = len(r["nodes"])
    for j, n in enumerate(r["nodes"]):
        VI[i, j] = v_of[n]
save["lig_vidx"] = VI
save["lig_nv"] = NV
save["n_v"] = np.array(len(allv))

np.savez(OUT / "ligaments.npz", **save)
print(f"\n写出 {OUT / 'ligaments.npz'}")
print(f"  韧带 {len(records)} 条，韧带节点 {len(allv)} 个，杆 {int(NTR.sum())} 根")

meta = [{k: r[k] for k in ("pid", "name", "ea", "thick", "a_bone", "b_bone",
                           "area", "L", "A_cross", "A_each")}
        | {"n_truss": len(r["trusses"]), "n_nodes": len(r["nodes"])}
        for r in records]
(OUT / "ligaments.json").write_text(
    json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
print(f"  元数据 -> ligaments.json")
