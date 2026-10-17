"""只读：皮质(hex) 与松质(tet) 是否共节点？—— 决定能否独立加密。"""
from pathlib import Path
import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
d = np.load(ROOT / "temp/thums/thums_lowerlimb.npz", allow_pickle=True)

BONES = {
    "tibia":   {"cort": ["81000700_solid"],
                "spon": ["81000600_solid", "81000601_solid"]},
    "fibula":  {"cort": ["81000900_solid"],
                "spon": ["81000800_solid", "81000801_solid"]},
    "talus":   {"cort": ["81001100_solid"],
                "spon": ["81001000_solid"]},
    "calc":    {"cort": ["81001300_solid"],
                "spon": ["81001200_solid"]},
}
KEYS = ("81000700_solid", "81000600_solid", "81000601_solid",
        "81000900_solid", "81000800_solid", "81000801_solid",
        "81001100_solid", "81001000_solid",
        "81001300_solid", "81001200_solid",
        "81001400_solid", "81001500_solid")


def unodes(keys):
    s = set()
    for k in keys:
        for row in d[k]:
            s.update(int(x) for x in row if x > 0)
    return s


print("=== 各 part 唯一节点数 ===")
for k in KEYS:
    print(f"  {k:20s} {len(unodes([k])):>6,}")

print("\n=== 皮质 vs 松质 节点重合 ===")
for b, g in BONES.items():
    c = unodes(g["cort"]); s = unodes(g["spon"])
    inter = c & s
    print(f"  {b:8s} cort {len(c):>5,} / spon {len(s):>6,} / 重合 {len(inter):>5,}"
          f"  ({len(inter)/max(len(c),1)*100:5.1f}% of cort, "
          f"{len(inter)/max(len(s),1)*100:5.1f}% of spon)")

print("\n=== 踝部 12 个 solid part 之间两两重合（看整体连通性）===")
sets = {k: unodes([k]) for k in KEYS}
names = list(KEYS)
hdr = "               " + " ".join(f"{n[2:7]:>5s}" for n in names)
print(hdr)
for i, a in enumerate(names):
    row = []
    for j, b in enumerate(names):
        ov = len(sets[a] & sets[b]) if i != j else len(sets[a])
        row.append(f"{ov:>5d}")
    print(f"  {a:>12s} " + " ".join(row))

# 全模型的节点复用情况：每个节点被几个 part 用
print("\n=== 节点复用统计 ===")
use = {}
for k in KEYS:
    for row in d[k]:
        for x in row:
            x = int(x)
            if x > 0:
                use.setdefault(x, set()).add(k)
from collections import Counter
cnt = Counter(len(v) for v in use.values())
print(f"  踝部用到节点 {len(use):,} 个；被 N 个 part 共用的分布: {dict(sorted(cnt.items()))}")
print(f"  全局节点 {d['node_ids'].size:,}；踝部未用到的 {d['node_ids'].size - len(use):,}")
