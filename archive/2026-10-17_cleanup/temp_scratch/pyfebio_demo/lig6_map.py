"""L2: 重新判定韧带附着骨（这次有舟骨了）→ 定出胫跟/胫舟/跟腓/下胫腓前。

同时做回归：确认原有 26 个部件与备份逐位一致。
"""
import numpy as np
from pathlib import Path

T = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
new = np.load(T / "thums_lowerlimb.npz", allow_pickle=True)
old = np.load(T / "thums_lowerlimb.npz.bak-preLIG", allow_pickle=True)

print("=== 回归：原 26 个部件是否与备份逐位一致 ===")
same = True
for k in old.files:
    if k in ("node_ids", "node_xyz"):
        continue
    a, b = old[k], new[k]
    ok = a.shape == b.shape and np.array_equal(a, b)
    if not ok:
        same = False
        print(f"  ❌ {k}: {a.shape} -> {b.shape}")
print(f"  {'✅ 全部一致' if same else '⚠️ 有差异'}")

nid, xyz = new["node_ids"], new["node_xyz"]
i_of = {int(v): i for i, v in enumerate(nid)}
print(f"\n节点 {len(nid):,}（备份 {len(old['node_ids']):,}）")

BONE = {
    "TIBIA": [81000700, 81000600], "FIBULA": [81000900, 81000800],
    "TALUS": [81001100, 81001000], "CALCANEUS": [81001300, 81001200],
    "NAVICULAR": [81001400, 81001401], "CUBOID": [81001500, 81001501],
    "MED_CUNEIF": [81001600, 81001601], "INT_CUNEIF": [81001700, 81001701],
    "LAT_CUNEIF": [81001800, 81001801],
}
bone_pts = {}
for nm, pids in BONE.items():
    pts = []
    for pid in pids:
        for kind in ("solid", "shell"):
            k = f"{pid:08d}_{kind}"
            if k in new:
                pts.append(xyz[[i_of[int(n)] for n in new[k].ravel()]])
    if pts:
        bone_pts[nm] = np.vstack(pts)

LIG = [
    (81100701, "TALOFIBULARE_POSTERIUS", "shell"),
    (81100801, "CALCANEOFIBULARE", "shell"),
    (81100901, "TALOFIBULARE_ANTERIUS", "shell"),
    (81101001, "TALOCALCANEUM_POSTERIUS", "shell"),
    (81101101, "TIBIOFIBULARE", "shell"),
    (81101201, "MEDIALE1", "shell"), (81101202, "MEDIALE2", "shell"),
    (81101203, "MEDIALE3", "shell"), (81101204, "MEDIALE4", "shell"),
    (81101401, "LONG_PLANTAR", "shell"),
    (81100100, "COLLAT_FIBULAR_TISSUES", "solid"),
    (81100200, "COLLAT_TIBIAL_TISSUES", "solid"),
] + [(81101301 + i, f"FOOT_{i+1}", "shell") for i in range(8)] + [
    (81101408, "FOOT_10", "shell"), (81101409, "FOOT_11", "shell"),
    (81101410, "CAPITIS_FIBULAE", "shell"),
]

print("\n=== 每条韧带：两端各归属哪块骨（阈值 2.0 mm）===")
RES = {}
for pid, nm, kind in LIG:
    k = f"{pid:08d}_{kind}"
    if k not in new:
        continue
    e = new[k]
    nodes = sorted({int(n) for n in e.ravel()})
    pts = xyz[[i_of[n] for n in nodes]]
    # 每个节点到各骨的最近距
    D = {b: np.min(np.linalg.norm(
        bone_pts[b][None, :, :] - pts[:, None, :], axis=2), axis=1)
        for b in bone_pts}
    # 归属：取最近且 < 2.0mm 的骨
    near = {}
    for b, d in D.items():
        m = d < 2.0
        if m.any():
            near[b] = int(m.sum())
    top = sorted(near.items(), key=lambda kv: -kv[1])
    RES[pid] = top
    span = np.linalg.norm(pts.max(0) - pts.min(0))
    tstr = ", ".join(f"{b}:{n}" for b, n in top if n >= 2) or "(无 <2mm 附着)"
    print(f"  {pid} {nm:26s} {len(e):3d}e {len(nodes):3d}n 跨度{span:5.1f}  {tstr}")

print("\n\n=== 与吴恺排序的对应（关键判断）===")
print("  吴恺: 胫跟 > 胫舟 > 跟腓 > 下胫腓前\n")
for pid, nm, _ in LIG:
    if pid not in RES:
        continue
    bs = {b for b, n in RES[pid]}
    if "CALCANEUS" in bs and "TIBIA" in bs:
        print(f"  ★ {pid} {nm}: TIBIA+CALCANEUS -> **胫跟 (tibiocalcaneal)**")
    if "NAVICULAR" in bs and "TIBIA" in bs:
        print(f"  ★ {pid} {nm}: TIBIA+NAVICULAR -> **胫舟 (tibionavicular)**")
    if "CALCANEUS" in bs and "FIBULA" in bs:
        print(f"  ★ {pid} {nm}: FIBULA+CALCANEUS -> **跟腓 (calcaneofibular)**")
    if "TIBIA" in bs and "FIBULA" in bs:
        print(f"  ★ {pid} {nm}: TIBIA+FIBULA -> **下胫腓 (syndesmosis)**")
    if "NAVICULAR" in bs and "CALCANEUS" in bs:
        print(f"  ★ {pid} {nm}: CALCANEUS+NAVICULAR -> **弹簧韧带 (spring)**")
