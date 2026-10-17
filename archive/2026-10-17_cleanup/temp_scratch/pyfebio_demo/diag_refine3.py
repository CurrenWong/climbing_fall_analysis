"""定位：hex8 加密后 4 个子单元里，究竟哪一个翻负、什么条件下翻。"""
import numpy as np
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis\temp\thums\cartilage_mesh.npz")
d = np.load(P, allow_pickle=True)
bpos = {int(v): d["bone_xyz"][i] for i, v in enumerate(d["bone_node_ids"])}
FAN8 = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
        (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]


def vol8(P8):
    return sum(float(np.dot(np.cross(P8[q[1]] - P8[q[0]], P8[q[2]] - P8[q[0]]),
                            P8[q[3]] - P8[q[0]])) / 6.0 for q in FAN8)


E = np.asarray(d["tibiotalar_tibia_elems"], dtype=np.int64)
off = np.asarray(d["tibiotalar_tibia_offsets"], float)
rows = [[int(x) for x in r if int(x) != 0] for r in E]
rows = [r for r in rows if len(r) == 8]

from collections import Counter
which = Counter()
tot = Counter()
for v in rows:
    P8 = np.array([bpos[x] if x > 0 else off[-x - 1] for x in v])
    b = [P8[i] for i in range(4)]
    o = [P8[4 + i] for i in range(4)]
    m = [(b[i] + b[(i + 1) % 4]) / 2 for i in range(4)]
    mo = [(o[i] + o[(i + 1) % 4]) / 2 for i in range(4)]
    cb, co = np.mean(b, 0), np.mean(o, 0)
    subs = [[b[0], m[0], cb, m[3], o[0], mo[0], co, mo[3]],
            [m[0], b[1], m[1], cb, mo[0], o[1], mo[1], co],
            [cb, m[1], b[2], m[2], co, mo[1], o[2], mo[2]],
            [m[3], cb, m[2], b[3], mo[3], co, mo[2], o[3]]]
    for k, s in enumerate(subs):
        tot[k] += 1
        if vol8(np.array(s)) <= 0:
            which[k] += 1

print("每个子单元位置的翻负计数:", dict(which), " / 总数:", dict(tot))

# 打印一个翻负例子，核对 base/outer 四边形的"同侧性"
for v in rows:
    P8 = np.array([bpos[x] if x > 0 else off[-x - 1] for x in v])
    b = [P8[i] for i in range(4)]
    o = [P8[4 + i] for i in range(4)]
    nb = np.cross(b[1] - b[0], b[2] - b[0])
    no = np.cross(o[1] - o[0], o[2] - o[0])
    mid = (np.mean(b, 0) + np.mean(o, 0)) / 2
    db = np.dot(nb, mid - np.mean(b, 0))
    do = np.dot(no, mid - np.mean(o, 0))
    if db * do > 0:            # 两个法向同侧 = base/outer 绕向"反了"
        print("\n★ 找到 base/outer 绕向不一致的单元")
        print(f"  base 法向·(中心-底心) = {db:+.3f}   outer = {do:+.3f}"
              f"   -> {'同侧(错)' if db*do>0 else '异侧(对)'}")
        for i in range(4):
            print(f"   b{i} {np.round(b[i],3)}   o{i} {np.round(o[i],3)}"
                  f"   d={np.round(o[i]-b[i],3)}")
        break

n_bad = 0
for v in rows:
    P8 = np.array([bpos[x] if x > 0 else off[-x - 1] for x in v])
    b = [P8[i] for i in range(4)]
    o = [P8[4 + i] for i in range(4)]
    mid = (np.mean(b, 0) + np.mean(o, 0)) / 2
    db = np.dot(np.cross(b[1] - b[0], b[2] - b[0]), mid - np.mean(b, 0))
    do = np.dot(np.cross(o[1] - o[0], o[2] - o[0]), mid - np.mean(o, 0))
    if db * do > 0:
        n_bad += 1
print(f"\nbase/outer 绕向不一致的 hex8 单元: {n_bad}/{len(rows)}")
