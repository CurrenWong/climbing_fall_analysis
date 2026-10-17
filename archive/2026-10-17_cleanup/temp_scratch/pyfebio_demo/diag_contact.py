"""查接触为何未检测：① FEBio 日志中的 contact 段 ② 两个接触面的法向朝向"""
import re
from pathlib import Path
import numpy as np

P = Path(r"D:\Project\climbing_fall_analysis")
log = (P / "temp/thums/febio/ankle_tibiotalar.log").read_text(encoding="latin-1",
                                                              errors="replace")
print("=== 日志中与 contact 相关的行 ===")
for l in log.split("\n"):
    if re.search(r"contact|interface|sliding|traction pair|penetrat", l, re.I):
        print("  " + l.strip())

# 复算两个接触面的定向
db = np.load(P / "temp/thums/thums_lowerlimb.npz")
nid, xyz = db["node_ids"], db["node_xyz"]
idx = {int(v): i for i, v in enumerate(nid)}
cg = np.load(P / "temp/thums/cartilage_mesh.npz")

# 骨节点 → 紧凑编号（与 build_feb 一致）
bone_ids = sorted({int(x) for j in ("tibiotalar_tibia", "tibiotalar_talus")
                   for row in cg[f"{j}_elems"] for x in row if int(x) > 0})
coords = [tuple(xyz[idx[v]]) for v in bone_ids]
off_of = {}
for j in ("tibiotalar_tibia", "tibiotalar_talus"):
    for k, p in enumerate(cg[f"{j}_offsets"]):
        off_of[f"off:{j}:{k}"] = len(coords)
        coords.append(tuple(p))
coords = np.asarray(coords, float)

print("\n=== 两个接触面（软骨外层）的法向统计 ===")
for j in ("tibiotalar_tibia", "tibiotalar_talus"):
    ns = []
    for row in cg[f"{j}_elems"]:
        vals = [int(x) for x in row if int(x) != 0]
        if not vals or vals[0] <= 0:
            continue
        m = len(vals) // 2
        pts = coords[[off_of[f"off:{j}:{-x - 1}"] for x in vals[m:]]]
        n = np.zeros(3)
        for i in range(len(pts)):
            n += np.cross(pts[i], pts[(i + 1) % len(pts)])
        ln = np.linalg.norm(n)
        if ln > 1e-12:
            ns.append(n / ln)
    ns = np.array(ns)
    mean = ns.mean(0)
    # 该关节另一侧骨的质心方向
    other = "81001000_solid" if j.endswith("tibia") else "81000600_solid"
    ob = db[other]
    oc = np.array([xyz[idx[int(v)]] for r in ob[:200] for v in r]).mean(0)
    # 本侧软骨 offset 点的质心
    own = np.array([coords[off_of[f"off:{j}:{k}"]] for k in range(len(cg[f"{j}_offsets"]))])
    ownc = own.mean(0)
    to_other = oc - ownc
    to_other = to_other / np.linalg.norm(to_other)
    print(f"  {j}: 面数 {len(ns)}")
    print(f"    平均法向       {np.round(mean, 4)}")
    print(f"    指向对侧骨方向 {np.round(to_other, 4)}")
    print(f"    点积 = {np.dot(mean, to_other):+.4f}  "
          f"({'朝外(正确)' if np.dot(mean, to_other) > 0 else '朝内(疑似反向!)'})")
