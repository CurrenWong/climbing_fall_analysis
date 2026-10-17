"""软骨 vs 骨：节点是否共享？几何是否重合？"""
import re
from pathlib import Path

import numpy as np

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio")


def load(fn):
    txt = Path(fn).read_text(encoding="latin-1")
    parts = {}
    for m in re.finditer(
            r'<Elements\b[^>]*type="(tet4|hex8|penta6)"[^>]*name="([^"]+)"[^>]*>(.*?)</Elements>',
            txt, re.S):
        nm = m.group(2)
        parts[nm] = {tuple(int(y) for y in s.split(","))
                     for s in re.findall(r"<elem\b[^>]*>\s*([\d,]+)\s*</elem>", m.group(3))}
    XYZ = {}
    for m in re.finditer(
            r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)', txt):
        XYZ[int(m.group(1))] = (float(m.group(2)), float(m.group(3)), float(m.group(4)))
    # NodeSet
    nss = {m.group(1): {int(x) for x in m.group(2).split(",") if x.strip()}
           for m in re.finditer(r'<NodeSet name="([^"]+)">([^<]*)</NodeSet>', txt)}
    # Surface
    sfs = {m.group(1): {int(y) for s in re.findall(r"<[a-z0-9]+ id=[^>]*>\s*([\d,]+)\s*</", m.group(2))
                        for y in s.split(",")}
           for m in re.finditer(r'<Surface name="([^"]+)"[^>]*>(.*?)</Surface>', txt, re.S)}
    return parts, XYZ, nss, sfs


for tag in ("_base.feb", "ankle_tibiotalar.feb"):
    p, XYZ, nss, sfs = load(W / tag)
    print("\n" + "=" * 74)
    print(f"{tag}")
    print("=" * 74)
    for nm in sorted(p):
        print(f"  {nm:26s} 单元 {len(p[nm]):>6d}  节点 {len({x for e in p[nm] for x in e}):>6d}")

    from scipy.spatial import cKDTree
    for side, bones, cart in (("胫骨", ("tibia_cort__hex8", "tibia_cort__tet4",
                                      "tibia_end_spon__tet4"),
                               "tibiotalar_tibia__tet4"),
                              ("距骨", ("talus_cort__hex8", "talus_spon__tet4"),
                               "tibiotalar_talus__tet4")):
        if cart not in p:
            continue
        bn = {x for b in bones if b in p for e in p[b] for x in e}
        cn = {x for e in p[cart] for x in e}
        inter = bn & cn
        B = np.array([XYZ[i] for i in sorted(bn)])
        C = np.array([XYZ[i] for i in sorted(cn)])
        d, _ = cKDTree(B).query(C)
        print(f"\n  [{side}] 骨节点 {len(bn)} / 软骨节点 {len(cn)}  "
              f"⇒ **共享 id {len(inter)}**")
        if inter:
            print(f"    共享示例 {sorted(inter)[:6]}")
        print(f"    软骨节点→最近骨节点距离：中位 {np.median(d):.4f} "
              f"p95 {np.percentile(d,95):.4f} max {d.max():.4f} mm"
              f"   重合(<1e-6) {int((d<1e-6).sum())}/{len(C)}")

    print("\n  NodeSet:")
    for k, v in nss.items():
        print(f"    {k:26s} {len(v)} 节点")
    print("  Surface:")
    for k, v in sfs.items():
        print(f"    {k:26s} {len(v)} 节点")
