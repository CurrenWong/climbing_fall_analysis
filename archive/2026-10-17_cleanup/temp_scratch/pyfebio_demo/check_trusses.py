"""检查生成的杆是否有零长/退化/重复（FEBio 遇零长杆会崩溃）。"""
import numpy as np
from pathlib import Path

OUT = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
db = np.load(OUT / "thums_lowerlimb.npz", allow_pickle=True)
lg = np.load(OUT / "ligaments.npz", allow_pickle=True)
nid, xyz = db["node_ids"], db["node_xyz"]
i_of = {int(v): i for i, v in enumerate(nid)}

n_lig = int(lg["n_lig"])
print(f"{'韧带':32s} {'杆数':>4} {'最短':>8} {'最长':>8} {'零长':>5} {'重复':>5}")
print("-" * 78)
bad = 0
allpairs = []
for i in range(n_lig):
    nm = str(lg["names"][i])
    tr = [[int(x) for x in t] for t in lg["truss"][i]
          if int(t[0]) > 0 and int(t[1]) > 0]
    L = []
    nzero = 0
    for a, b in tr:
        d = float(np.linalg.norm(xyz[i_of[a]] - xyz[i_of[b]]))
        L.append(d)
        if d < 1e-6:
            nzero += 1
    dup = len(tr) - len({tuple(t) for t in tr})
    bad += nzero
    allpairs += [(nm, a, b) for a, b in tr]
    print(f"{nm:32s} {len(tr):>4} {min(L) if L else 0:>8.3f} "
          f"{max(L) if L else 0:>8.3f} {nzero:>5} {dup:>5}")

print(f"\n零长杆合计 {bad}")
print(f"全部杆 {len(allpairs)}")

# 重复节点检查
from collections import Counter
c = Counter()
for nm, a, b in allpairs:
    c[(a, b)] += 1
rep = {k: v for k, v in c.items() if v > 1}
print(f"跨韧带重复的 (a,b) 对: {len(rep)}")
if rep:
    for (a, b), v in list(rep.items())[:8]:
        print(f"   {a},{b} x{v}")

# 每根杆的两端是否真在不同骨上
print("\n=== 检查是否有杆的两端落在同一骨 ===")
BONE = {}
for pid, dom in ((81000700, "tibia"), (81000600, "tibia"),
                 (81000601, "tibia"), (81000900, "fibula"),
                 (81000800, "fibula"), (81001100, "talus"),
                 (81001000, "talus"), (81001300, "calcaneus"),
                 (81001200, "calcaneus"), (81001400, "navicular")):
    for kind in ("solid", "shell"):
        k = f"{pid:08d}_{kind}"
        if k in db:
            BONE.setdefault(dom, set()).update({int(v) for v in db[k].ravel()})

owner = {}
for d, s in BONE.items():
    for n in s:
        owner.setdefault(n, []).append(d)

same = []
for nm, a, b in allpairs:
    oa, ob = set(owner.get(a, [])), set(owner.get(b, []))
    if oa & ob:
        same.append((nm, a, b, sorted(oa & ob)))
print(f"两端同骨的杆: {len(same)}")
for s in same[:10]:
    print("   ", s)
