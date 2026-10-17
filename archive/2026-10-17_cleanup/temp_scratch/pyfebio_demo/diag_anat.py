"""软骨结构解剖：哪个面共享骨节点？哪个面朝关节？漂移从哪来？"""
import re
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio")
FN = W / "_base.feb"
txt = FN.read_text(encoding="latin-1")

parts = {}
for m in re.finditer(
        r'<Elements\b[^>]*type="(tet4|hex8|penta6)"[^>]*name="([^"]+)"[^>]*>(.*?)</Elements>',
        txt, re.S):
    parts[m.group(2)] = {int(y) for s in
                         re.findall(r"<elem\b[^>]*>\s*([\d,]+)\s*</elem>", m.group(3))
                         for y in s.split(",")}
XYZ = {}
for m in re.finditer(
        r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)', txt):
    XYZ[int(m.group(1))] = np.array([float(m.group(2)), float(m.group(3)),
                                     float(m.group(4))])
sfs = {}
for m in re.finditer(r'<Surface name="([^"]+)"[^>]*>(.*?)</Surface>', txt, re.S):
    sfs[m.group(1)] = {int(y) for s in
                       re.findall(r"<[a-z0-9]+ id=[^>]*>\s*([\d,]+)\s*</", m.group(2))
                       for y in s.split(",")}

AX = np.array([0.8153, 0.135, 0.5631])
AX = AX / np.linalg.norm(AX)


def pos(i):
    return float(XYZ[i] @ AX)


for side, bone, cart in (("胫骨", ("tibia_cort__hex8", "tibia_cort__tet4",
                                 "tibia_end_spon__tet4", "tibia_marrow__tet4"),
                          "tibiotalar_tibia__tet4"),
                         ("距骨", ("talus_cort__hex8", "talus_spon__tet4"),
                          "tibiotalar_talus__tet4")):
    bn = {x for b in bone if b in parts for x in parts[b]}
    cn = parts[cart]
    s_outer = sfs[f"surf_{cart.replace('__tet4','')}"]
    s_boneface = sfs[f"talus_cort_bottom"] if side == "距骨" else None
    print("=" * 76)
    print(f"【{side}】骨 {len(bn)} 节点 / 软骨 {len(cn)} 节点")
    print("=" * 76)
    shared = bn & cn
    print(f"  软骨节点中属于骨的（共享 id）: {len(shared)}")
    print(f"  surface surf_{cart.replace('__tet4','')}: {len(s_outer)} 节点")
    print(f"    ⊂ 骨? {s_outer <= bn}   ⊂ 软骨? {s_outer <= cn}   "
          f"∩骨={len(s_outer & bn)}")
    print(f"  软骨节点 − 该 surface = {len(cn - s_outer)} 个（另一侧）")

    # 沿轴位置
    pb = np.array([pos(i) for i in bn])
    po = np.array([pos(i) for i in s_outer])
    pc = np.array([pos(i) for i in cn])
    pr = np.array([pos(i) for i in (cn - s_outer)])

    def rng(a):
        return f"[{a.min():8.3f},{a.max():8.3f}] 中位{a.mean():8.3f}"
    print(f"\n  沿加载轴位置（越大越靠胫骨近端）:")
    print(f"    骨节点         {rng(pb)}")
    print(f"    软骨全部       {rng(pc)}")
    print(f"    surf_*(接触面) {rng(po)}   ← {'靠骨侧' if po.mean() > pc.mean() else '靠关节侧'}")
    print(f"    软骨另一侧     {rng(pr)}   ← {'靠关节侧' if pr.mean() < po.mean() else '靠骨侧'}")

    # 语义确认：软骨厚度方向
    print(f"\n  软骨厚度（沿轴） ≈ {abs(po.mean()-pr.mean()):.3f} mm；"
          f"软骨到骨面间隙 = {pc.min()-pb.max():+.3f} mm" if side == "胫骨" else "")
