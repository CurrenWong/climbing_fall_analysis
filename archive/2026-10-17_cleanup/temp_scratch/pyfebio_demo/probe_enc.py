"""查清 {j}_elems 的真实编码：正负号分布、底面是否真是骨节点、共享度。"""
import numpy as np
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis\temp\thums\cartilage_mesh.npz")
d = np.load(P, allow_pickle=True)
bone_ids = set(int(v) for v in d["bone_node_ids"])
n_off = {j: len(d[f"{j}_offsets"]) for j in
         ("tibiotalar_tibia", "tibiotalar_talus",
          "subtalar_talus", "subtalar_calcaneus")}

for j, nn in n_off.items():
    E = np.asarray(d[f"{j}_elems"], dtype=np.int64)
    print(f"\n=== {j}  offsets={nn}  elems={E.shape} ===")
    widths = {}
    signpat = {}
    for row in E:
        v = [int(x) for x in row if int(x) != 0]
        w = len(v)
        widths[w] = widths.get(w, 0) + 1
        pat = "".join("+" if x > 0 else "-" for x in v)
        signpat[pat] = signpat.get(pat, 0) + 1
    print("  宽度分布:", widths)
    print("  符号模式:", dict(sorted(signpat.items(), key=lambda kv: -kv[1])))

    # 逐宽度看正/负引用是否越界
    for w in widths:
        rows = [[int(x) for x in r if int(x) != 0] for r in E]
        rows = [r for r in rows if len(r) == w]
        pos = [x for r in rows for x in r if x > 0]
        neg = [x for r in rows for x in r if x < 0]
        pb = sum(1 for x in pos if x in bone_ids)
        ob = max((-x - 1) for x in neg) if neg else -1
        print(f"   w={w}: n={len(rows)}  正引用 {len(pos)} 个, 其中属于骨ID {pb} 个"
              f" | 负引用 {len(neg)} 个, 最大下标 {ob} (应 < {nn})")

    # 底面(前一半)的节点共享度
    for w in (8, 6):
        rows = [[int(x) for x in r if int(x) != 0] for r in E]
        rows = [r for r in rows if len(r) == w]
        if not rows:
            continue
        h = w // 2
        base_u = {x for r in rows for x in r[:h]}
        out_u = {x for r in rows for x in r[h:]}
        print(f"   w={w} 底面唯一节点 {len(base_u)}（{h}*{len(rows)}={h*len(rows)} 引用）"
              f" | 外表面唯一节点 {len(out_u)}")
        print(f"        底面里属于骨ID的: "
              f"{sum(1 for x in base_u if x > 0 and x in bone_ids)}"
              f" / 负引用: {sum(1 for x in base_u if x < 0)}")
