"""①关节面初始间隙/穿透  ②高位移骨节点是否就是关节界面节点。"""
import inspect
import re
import textwrap
from pathlib import Path

import h5py
import numpy as np
import pyfebio.xplt as X
from scipy.spatial import cKDTree

ROOT = Path(r"D:\Project\climbing_fall_analysis")
W = ROOT / "temp/thums/febio"
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


def nodes_of(dom):
    return parts[dom]


P = np.array([XYZ[i] for i in sorted(sfs["surf_tibiotalar_tibia"])])
Q = np.array([XYZ[i] for i in sorted(sfs["surf_tibiotalar_talus"])])

print("=" * 76)
print("① 关节两软骨面 初始最小距离（穿透 ⇒ 惩罚力爆炸）")
print("=" * 76)
d, idx = cKDTree(Q).query(P)
print(f"  胫骨软骨关节面 → 距骨软骨关节面：min {d.min():.4f}  "
      f"p5 {np.percentile(d,5):.4f} 中位 {np.median(d):.4f} max {d.max():.4f} mm")
print(f"  'min<0.1mm' 的面节点数 = {int((d<0.1).sum())}/{len(P)}")
print(f"  >0.5mm 的 = {int((d>0.5).sum())}   >1mm 的 = {int((d>1.0).sum())}")
if (d < 0.5).any():
    i = int(np.argmin(d))
    print(f"  最近一对: 胫骨 node {sorted(sfs['surf_tibiotalar_tibia'])[i]} "
          f"{P[i]} ↔ 距骨 {Q[idx[i]]}  d={d[i]:.4f}")

# bbox
print("\n  两关节面 bbox（沿加载轴）：")
AX = np.array([0.8153, 0.135, 0.5631])
AX = AX / np.linalg.norm(AX)
print(f"    胫骨软骨面 轴投影 [{P@AX:.0f}]".replace(".0", "", 1) if False else
      f"    胫骨软骨面 轴投影 [{ (P@AX).min():.2f}, {(P@AX).max():.2f}]")
print(f"    距骨软骨面 轴投影 [{ (Q@AX).min():.2f}, {(Q@AX).max():.2f}]")
print(f"    ⇒ 重叠深度 = {min((P@AX).max(), (Q@AX).max()) - max((P@AX).min(), (Q@AX).min()):.3f} mm")

# ── ② 高位移骨节点 ──
for _o, _n in (('domain["name"]', 'domain.get("name", str(domain["id"][0]))'),
               ('surface["name"]', 'surface.get("name", str(surface["id"][0]))'),
               ("{domain['name']}", "{domain.get('name', str(domain['id'][0]))}"),
               ("{surface['name']}", "{surface.get('name', str(surface['id'][0]))}"),
               ('if nodeset["name"] == "":', 'if nodeset.get("name", "") == "":'),
               ("{nodeset['name']}", "{nodeset.get('name', str(nodeset['id'][0]))}"),
               ("{elementset['name']}",
                "{elementset.get('name', str(elementset['id'][0]))}")):
    pass
_src = textwrap.dedent(inspect.getsource(X.parse_mesh))
for _o, _n in (('domain["name"]', 'domain.get("name", str(domain["id"][0]))'),
               ('surface["name"]', 'surface.get("name", str(surface["id"][0]))'),
               ("{domain['name']}", "{domain.get('name', str(domain['id'][0]))}"),
               ("{surface['name']}", "{surface.get('name', str(surface['id'][0]))}"),
               ('if nodeset["name"] == "":', 'if nodeset.get("name", "") == "":'),
               ("{nodeset['name']}", "{nodeset.get('name', str(nodeset['id'][0]))}"),
               ("{elementset['name']}",
                "{elementset.get('name', str(elementset['id'][0]))}")):
    _src = _src.replace(_o, _n)
_ns = X.__dict__
exec(compile(_src, "<p>", "exec"), _ns)
X.parse_mesh = _ns["parse_mesh"]
_ps = textwrap.dedent(inspect.getsource(X.parse_state))
for _a, _b in (
    ('set_name = mesh_dict[set_lut[key]][set_id + 1]',
     'set_name = mesh_dict[set_lut[key]].get(set_id + 1) or '
     'mesh_dict[set_lut[key]].get(int(set_id) + 1) or f"set{set_id}"'),
    ('set_name = mesh_dict[set_lut[key]][set_id]',
     'set_name = mesh_dict[set_lut[key]].get(set_id) or '
     'mesh_dict[set_lut[key]].get(int(set_id)) or f"set{set_id}"'),
):
    _ps = _ps.replace(_a, _b)
_ns2 = X.__dict__
exec(compile(_ps, "<p2>", "exec"), _ns2)
X.parse_state = _ns2["parse_state"]
X.to_hdf5(str(W / "_base.xplt"), str(W / "_ba.h5"))
with h5py.File(W / "_ba.h5", "r") as h:
    n = len(h["states"]) - 1
    U = np.asarray(h[f"states/{n}/node_data/displacement/all"][()], float)
mag = np.linalg.norm(U, axis=1)

print("\n" + "=" * 76)
print("② 高位移节点归属（末帧 |u|>0.6mm）")
print("=" * 76)
cart_shared = (nodes_of("tibia_cort__hex8") | nodes_of("tibia_cort__tet4")
               | nodes_of("tibia_end_spon__tet4") | nodes_of("tibia_marrow__tet4")) \
    & nodes_of("tibiotalar_tibia__tet4")
hi = {i + 1 for i in range(len(mag)) if mag[i] > 0.6}
print(f"  |u|>0.6mm 的节点 {len(hi)} 个")
for dom, ns in sorted(parts.items(), key=lambda kv: -len(kv[1])):
    k = len(hi & ns)
    if k:
        print(f"    {dom:26s} {k:>5d} / {len(ns):>6d}")
print(f"  其中属于「胫骨↔软骨共享界面」的: {len(hi & cart_shared)} / "
      f"{len(hi)}（共享界面共 {len(cart_shared)} 节点）")
# 胫骨骨内部分布
tb = nodes_of("tibia_cort__hex8")
tb_mag = np.array([mag[i - 1] for i in sorted(tb)])
print(f"\n  胫骨皮质 |u| 分位: p50={np.percentile(tb_mag,50):.4f} "
      f"p90={np.percentile(tb_mag,90):.4f} p99={np.percentile(tb_mag,99):.4f} "
      f"max={tb_mag.max():.4f}")
# 高位移节点的轴位置（是只在关节端还是遍布整根骨）
tbn = sorted(tb)
for lo, hi_ in ((0, 0.3), (0.3, 0.6), (0.6, 0.9), (0.9, 3.0)):
    sel = [i for i in tbn if lo <= mag[i - 1] < hi_]
    if sel:
        pr = np.array([XYZ[i] @ AX for i in sel])
        print(f"    |u|∈[{lo},{hi_}) : {len(sel):>5d} 节点  轴位置中位 {np.median(pr):.1f}  "
              f"范围[{pr.min():.1f},{pr.max():.1f}]")
