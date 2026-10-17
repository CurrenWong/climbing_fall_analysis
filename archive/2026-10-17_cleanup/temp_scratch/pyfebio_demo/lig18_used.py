"""精确重建 build_feb.py 的 `used` 集合，检查韧带附着节点是否已在网格中。"""
import numpy as np
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
OUT = ROOT / "temp/thums"
db = np.load(OUT / "thums_lowerlimb.npz", allow_pickle=True)
nid, xyz = db["node_ids"], db["node_xyz"]

BONE_PARTS = {
    "tibia_cort": "81000700_solid",
    "tibia_end_spon": "81000600_solid",
    "tibia_marrow": "81000601_solid",
    "talus_cort": "81001100_solid",
    "talus_spon": "81001000_solid",
}
# FULL_ANKLE 模式还含跟骨
for k, v in (("calcaneus_cort", "81001300_solid"),
             ("calcaneus_spon", "81001200_solid")):
    BONE_PARTS[k] = v

used = set()
for dom, key in BONE_PARTS.items():
    for row in db[key]:
        conn = []
        for x in row:
            xi = int(x)
            if xi not in (0, -1) and xi not in conn:
                conn.append(xi)
        used.update(conn)
print(f"重建 used = {len(used):,} 节点")

# 各骨部件自身的节点
part_nodes = {k: {int(v) for v in db[k].ravel()} for k in BONE_PARTS.values()}

LIG = [(81100701, "PTFL 距腓后", 180), (81100801, "CFL 跟腓", 510),
       (81100901, "ATFL 距腓前", 180), (81101001, "PTCL 距跟后", 75.4),
       (81101101, "下胫腓", 75.4), (81101201, "内侧1", 180),
       (81101202, "内侧2 胫跟", 510), (81101203, "内侧3 胫舟", 510),
       (81101204, "内侧4", 180)]

print(f"\n{'pid':>9} {'名称':14s} {'总节点':>5} {'在used':>7} {'附着':>5} "
      f"{'附着在used':>9} {'挂在哪个骨部件':>28}")
print("-" * 92)
for pid, nm, ea in LIG:
    nodes = {int(v) for v in db[f"{pid:08d}_shell"].ravel()}
    inu = nodes & used
    att = {}
    for k, pn in part_nodes.items():
        inter = nodes & pn
        if inter:
            att[k] = inter
    att_in_used = {k: len(v & used) for k, v in att.items()}
    att_all = sum(len(v) for v in att.values())
    att_all_used = sum(att_in_used.values())
    s = ", ".join(f"{k}:{n}/{len(att[k])}" for k, n in att_in_used.items()
                  if len(att[k]) >= 2)
    print(f"{pid:>9} {nm:14s} {len(nodes):>5} {len(inu):>7} {att_all:>5} "
          f"{att_all_used:>9}  {s}")
