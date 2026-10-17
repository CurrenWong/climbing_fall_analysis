"""决定性验证：R_LIG_FOOT 1..9 (81101301-309) 是否才是与 R_ 骨匹配的踝韧带？"""
import numpy as np
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
z = np.load(ROOT / "temp/thums/thums_lowerlimb.npz", allow_pickle=True)
print("npz 键（含 811013xx 的）:", [k for k in z.keys() if k.startswith("811013")][:40])
print("npz 键（含 811014xx 的）:", [k for k in z.keys() if k.startswith("811014")][:40])
print("npz 键数量:", len(z.keys()))

if "node_ids" not in z:
    raise SystemExit("无 node_ids")
nid, xyz = z["node_ids"], z["node_xyz"]
i_of = {int(v): i for i, v in enumerate(nid)}

BONES = {
    "tibia": [81000700, 81000600, 81000601],
    "fibula": [81000900, 81000800],
    "talus": [81001100, 81001000],
    "calcaneus": [81001300, 81001200],
    "navicular": [81001400, 81001401],
}
BONE_IDS = {}
for b, pids in BONES.items():
    s = set()
    for pid in pids:
        for kind in ("solid", "shell"):
            k = f"{pid:08d}_{kind}"
            if k in z:
                s |= {int(v) for v in z[k].ravel()}
    BONE_IDS[b] = s

CAND = [(81101301, "R_LIG_FOOT 1"), (81101302, "R_LIG_FOOT 2"),
        (81101303, "R_LIG_FOOT 3"), (81101304, "R_LIG_FOOT 4"),
        (81101305, "R_LIG_FOOT 5"), (81101306, "R_LIG_FOOT 6"),
        (81101307, "R_LIG_FOOT 7"), (81101308, "R_LIG_FOOT 8"),
        (81101309, "R_LIG_FOOT 9"),
        (81100901, "L_ATFL(现用)"), (81100801, "L_CFL(现用)"),
        (81100701, "L_PTFL(现用)")]

print(f"\n{'pid':>10} {'名':<16} {'壳元':>5} {'节点':>5}  各骨命中数(tibia/fibula/talus/calcn/navic)  "
      f"{'最大两骨距':>10}")
print("─" * 108)
for pid, nm in CAND:
    key = f"{pid:08d}_shell"
    if key not in z:
        print(f"{pid:>10} {nm:<16} 无 shell")
        continue
    e = z[key]
    nodes = sorted({int(v) for v in e.ravel()})
    cnt = {b: sum(1 for n in nodes if n in BONE_IDS[b]) for b in BONES}
    # 主要两骨（命中最多的两个）的质心距
    top = sorted(cnt.items(), key=lambda kv: -kv[1])[:2]
    P = xyz[[i_of[n] for n in nodes]]
    info = ""
    if all(v > 0 for _, v in top):
        (b1, _), (b2, _) = top
        c1 = xyz[[i_of[n] for n in nodes if n in BONE_IDS[b1]]].mean(0)
        c2 = xyz[[i_of[n] for n in nodes if n in BONE_IDS[b2]]].mean(0)
        info = f"{b1}↔{b2} {np.linalg.norm(c2-c1):.2f} mm"
    print(f"{pid:>10} {nm:<16} {len(e):>5} {len(nodes):>5}  "
          f"{cnt['tibia']:>5}/{cnt['fibula']:>5}/{cnt['talus']:>5}/"
          f"{cnt['calcaneus']:>5}/{cnt['navicular']:>5}   {info}")
