"""决定性检查：当前 .feb 网格里是否存在各韧带的附着节点（按坐标匹配）。"""
import numpy as np
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
FEB = ROOT / "temp/thums/febio/ankle_tibiotalar.feb"
z = np.load(ROOT / "temp/thums/thums_lowerlimb.npz", allow_pickle=True)
nid, xyz = z["node_ids"], z["node_xyz"]
i_of = {int(v): i for i, v in enumerate(nid)}

# 解析 .feb 的 <Nodes name="all">
coords = []
in_nodes = False
for line in FEB.read_text(encoding="utf-8").split("\n"):
    s = line.strip()
    if "<Nodes name=\"all\">" in s:
        in_nodes = True
        continue
    if in_nodes and s.startswith("</Nodes>"):
        break
    if in_nodes and s.startswith("<node"):
        v = s.split(">")[1].split("<")[0]
        coords.append([float(x) for x in v.split(",")])
C = np.asarray(coords)
print(f".feb 节点数 {len(C):,}")

# 坐标->feb序号 的查找（四舍五入到 1e-6 匹配）
key = {}
for i, p in enumerate(C):
    key[tuple(np.round(p, 5))] = i + 1

BONE_PIDS = [81000700, 81000600, 81000900, 81000800, 81001100, 81001000,
             81001300, 81001200, 81001400, 81001401]
bone_nodes = set()
for pid in BONE_PIDS:
    for kind in ("solid", "shell"):
        k = f"{pid:08d}_{kind}"
        if k in z:
            bone_nodes |= {int(v) for v in z[k].ravel()}

LIG = [(81100701, "PTFL 距腓后"), (81100801, "CFL 跟腓"),
       (81100901, "ATFL 距腓前"), (81101001, "PTCL 距跟后"),
       (81101101, "下胫腓"), (81101201, "内侧1"),
       (81101202, "内侧2 胫跟"), (81101203, "内侧3 胫舟"),
       (81101204, "内侧4"), (81101401, "跖长")]

print(f"\n{'pid':>9} {'名称':14s} {'总节点':>5} {'连骨':>5} {'连骨且在feb':>10} "
      f"{'缺节点数':>8} {'→ 新增?':>8}")
print("-" * 90)
tot_missing = 0
for pid, nm in LIG:
    e = z[f"{pid:08d}_shell"]
    nodes = {int(v) for v in e.ravel()}
    attached = nodes & bone_nodes
    in_feb = [n for n in attached
              if tuple(np.round(xyz[i_of[n]], 5)) in key]
    missing = [n for n in nodes
               if tuple(np.round(xyz[i_of[n]], 5)) not in key]
    tot_missing += len(missing)
    print(f"{pid:>9} {nm:14s} {len(nodes):>5} {len(attached):>5} "
          f"{len(in_feb):>10} {len(missing):>8} {'是' if missing else '否':>8}")

print(f"\n所有韧带节点中，不在 .feb 里的合计 {tot_missing} 个（需新增节点）")
print(f"当前 .feb 单元/节点：见 log；新增节点后需重算 unit 总数")
