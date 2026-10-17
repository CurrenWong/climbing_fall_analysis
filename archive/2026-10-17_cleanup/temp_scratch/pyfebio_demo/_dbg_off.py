import numpy as np

cg = np.load(r"D:\Project\climbing_fall_analysis\temp\thums\cartilage_mesh.npz")
for j in ("tibiotalar_tibia", "tibiotalar_talus",
          "subtalar_talus", "subtalar_calcaneus"):
    e, o = cg[f"{j}_elems"], cg[f"{j}_offsets"]
    lens = sorted({len([x for x in row if int(x) != 0]) for row in e})
    negs = [int(x) for row in e for x in row if int(x) < 0]
    print(f"{j:22s} elems={e.shape} offsets={len(o):4d} "
          f"vals_len={lens} neg_range={min(negs)}..{max(negs)} "
          f"-> off_idx 0..{-min(negs)-1}")
