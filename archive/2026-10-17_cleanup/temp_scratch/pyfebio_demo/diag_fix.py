"""固定集到底落在哪？胫骨↔距骨是否共享节点？（这是载荷路径的真相）"""
import re
from pathlib import Path

import numpy as np

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio")
AX = np.array([0.8153, 0.135, 0.5631])
AX = AX / np.linalg.norm(AX)

for tag in ("B", "C", "E"):
    p = W / f"_{tag}.feb"
    if not p.exists():
        continue
    txt = p.read_text(encoding="latin-1")
    els = {}
    for m in re.finditer(
            r'<Elements\b[^>]*type="(tet4|hex8|penta6)"[^>]*name="([^"]+)"[^>]*>(.*?)</Elements>',
            txt, re.S):
        els[m.group(2)] = {int(y) for s in
                           re.findall(r"<elem\b[^>]*>\s*([\d,]+)\s*</elem>", m.group(3))
                           for y in s.split(",")}
    nss = {m.group(1): {int(x) for x in m.group(2).split(",") if x.strip()}
           for m in re.finditer(r'<NodeSet name="([^"]+)">([^<]*)</NodeSet>', txt)}
    XYZ = {}
    for m in re.finditer(
            r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)', txt):
        XYZ[int(m.group(1))] = np.array([float(m.group(2)), float(m.group(3)),
                                         float(m.group(4))])
    print("\n" + "=" * 78)
    print(f"{tag}.feb")
    print("=" * 78)
    # 固定集的节点在哪些部件里
    fix = nss["talus_cort_bottom_ns"]
    print(f"  talus_cort_bottom_ns {len(fix)} 节点，落在：")
    for d, s in els.items():
        k = len(fix & s)
        if k:
            print(f"    {d:28s} {k:>5d}/{len(fix)}")
    pr = np.array([XYZ[i] @ AX for i in sorted(fix)])
    print(f"    轴位置 [{pr.min():.2f},{pr.max():.2f}] 中位 {np.median(pr):.2f}")
    # 各骨域节点数
    for d in ("tibia_cort__hex8", "talus_cort__hex8"):
        print(f"  {d:26s} {len(els.get(d,[]))} 节点")
    ti, ta = els.get("tibia_cort__hex8", set()), els.get("talus_cort__hex8", set())
    print(f"  ⚠ tibia_cort__hex8 ∩ talus_cort__hex8 = {len(ti & ta)} 节点")
    # 全部骨的并集 vs 全部软骨
    allbones = set().union(*[s for d, s in els.items() if "cart" not in d
                             and "__tet4" not in d or "cort" in d or "spon" in d
                             or "marrow" in d])
    allbones = set().union(*[s for d, s in els.items()
                             if not d.startswith("tibiotalar")])
    carts = set().union(*[s for d, s in els.items() if d.startswith("tibiotalar")])
    print(f"  全部骨节点 {len(allbones)}  全部软骨节点 {len(carts)}  "
          f"共享 {len(allbones & carts)}")
    # 软骨的「底面」节点（共享骨）在哪
    for c in ("tibiotalar_tibia__tet4", "tibiotalar_talus__tet4"):
        if c in els:
            sh = els[c] & allbones
            print(f"    {c:26s} {len(els[c])} 节点，其中 {len(sh)} 个是骨节点")
            if sh:
                q = np.array([XYZ[i] @ AX for i in sh])
                print(f"        共享节点轴位置 [{q.min():.2f},{q.max():.2f}] "
                      f"中位 {np.median(q):.2f}")
    # 两骨是否连通：从固定集出发看能到哪
    print(f"  fix 集也在 tibia 里？ {len(fix & els.get('tibia_cort__hex8', set()))} 个")
