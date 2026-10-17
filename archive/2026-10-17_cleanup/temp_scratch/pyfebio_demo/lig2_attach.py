"""L1/L2: 韧带几何 + 附着骨判定（每个韧带的两端各贴近哪块骨）。"""
import json
import numpy as np
from pathlib import Path

T = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
z = np.load(T / "thums_lowerlimb.npz", allow_pickle=True)
nid, xyz = z["node_ids"], z["node_xyz"]
i_of = {int(v): i for i, v in enumerate(nid)}

LIG = {
    81100701: "LIG_TALOFIBULARE_POSTERIUS",
    81100801: "LIG_CALCANEOFIBULARE(跟腓)",
    81100901: "LIG_TALOFIBULARE_ANTERIUS",
    81101001: "LIG_TALOCALCANEUM_POSTERIUS",
    81101101: "LIG_TIBIOFIBULARE(下胫腓)",
    81101201: "LIG_MEDIALE1",
    81101202: "LIG_MEDIALE2",
    81101203: "LIG_MEDIALE3",
    81101204: "LIG_MEDIALE4",
    81101401: "LIG_LONG_PLANTAR",
}
BONE = {
    "R_TIBIA_CORT": 81000700,
    "right_tibia_end_spon": 81000600,
    "R_FIBULA_CORT": 81000900,
    "right_fibula_end_spon": 81000800,
    "R_TALUS_CORT": 81001100,
    "R_TALUS_SPON": 81001000,
    "R_CALCANEUS_CORT": 81001300,
    "R_CALCANEUS_SPON": 81001200,
}

# 骨表面点集（皮质+松质全用）
bone_pts = {}
for nm, pid in BONE.items():
    k = f"{pid:08d}_solid"
    if k in z:
        bone_pts[nm] = xyz[[i_of[int(n)] for n in z[k].ravel()]]

print("=== 各韧带：节点数 / 单元数 / 长度 / 最贴近的两块骨 ===")
for pid, nm in LIG.items():
    e = z[f"{pid:08d}_shell"]
    nodes = sorted({int(n) for n in e.ravel()})
    pts = xyz[[i_of[n] for n in nodes]]
    # 每个韧带节点到各骨的最小距 -> 统计"最近骨"
    assign = {}
    for j, p in enumerate(pts):
        bn = min(bone_pts, key=lambda b: np.min(
            np.linalg.norm(bone_pts[b] - p, axis=1)))
        assign.setdefault(bn, []).append(j)
    # 韧带整体跨度
    span = np.linalg.norm(pts.max(0) - pts.min(0))
    top = sorted(assign.items(), key=lambda kv: -len(kv[1]))
    tstr = ", ".join(f"{b}x{len(v)}" for b, v in top[:3])
    print(f"\n  {pid} {nm}")
    print(f"    {len(e)} 单元 / {len(nodes)} 节点  包围盒对角 {span:.1f} mm")
    print(f"    节点归属: {tstr}")

print("\n\n=== summary.json 里有没有舟骨/其它跗骨？===")
s = json.loads((T / "summary.json").read_text(encoding="utf-8"))
print("keys:", list(s.keys())[:20])
for k in ("part_titles", "parts", "found", "targets"):
    if k in s:
        v = s[k]
        print(f"\n{s[k] if not isinstance(v, dict) else ''}")
        if isinstance(v, dict):
            hits = {kk: vv for kk, vv in v.items()
                    if any(t in str(vv).upper() for t in
                           ("NAVIC", "CUBOID", "CUNEIF", "TALUS", "CALCAN",
                            "LIG", "TIBIA", "FIBULA"))}
            print(f"  {len(v)} 条，含关键词 {len(hits)} 条:")
            for kk, vv in list(hits.items())[:40]:
                print(f"    {kk}: {vv}")
        break
