"""tie 是否真抓住：界面滑移量（变形后软骨底面到骨面的距离）+ 软骨体积对比。"""
import inspect
import re
import textwrap
from pathlib import Path

import h5py
import numpy as np
import pyfebio.xplt as X
from scipy.spatial import cKDTree

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio")
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


def parse(fn):
    txt = Path(fn).read_text(encoding="latin-1")
    els = {}
    for m in re.finditer(
            r'<Elements\b[^>]*type="(tet4|hex8|penta6)"[^>]*name="([^"]+)"[^>]*>(.*?)</Elements>',
            txt, re.S):
        els[m.group(2)] = ([tuple(int(y) for y in s.split(",")) for s in
                            re.findall(r"<elem\b[^>]*>\s*([\d,]+)\s*</elem>", m.group(3))],
                           m.group(1))
    XYZ = {}
    for m in re.finditer(
            r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)', txt):
        XYZ[int(m.group(1))] = np.array([float(m.group(2)), float(m.group(3)),
                                         float(m.group(4))])
    sf = {}
    for m in re.finditer(r'<Surface name="([^"]+)"[^>]*>(.*?)</Surface>', txt, re.S):
        sf[m.group(1)] = [tuple(int(y) for y in s.split(",")) for s in
                          re.findall(r"<[a-z0-9]+ id=[^>]*>\s*([\d,]+)\s*</", m.group(2))]
    return els, XYZ, sf


def tet_vol(p):
    a, b, c, d = p
    return abs(np.dot(b - a, np.cross(c - a, d - a))) / 6.0


print("=" * 76)
print("① 软骨体积（baseline vs m2）")
print("=" * 76)
vol = {}
for k in ("B", "E"):
    els, XYZ, sf = parse(W / f"_{k}.feb")
    for dom in ("tibiotalar_tibia__tet4", "tibiotalar_talus__tet4"):
        v = sum(tet_vol(np.array([XYZ[i] for i in e])) for e in els[dom][0])
        vol[(k, dom)] = (v, len(els[dom][0]))
print(f"  {'部件':26s} {'baseline(B)':>22} {'m2(E)':>22} {'倍数':>7}")
for dom in ("tibiotalar_tibia__tet4", "tibiotalar_talus__tet4"):
    vb, nb = vol[("B", dom)]
    ve, ne = vol[("E", dom)]
    print(f"  {dom:26s} {vb:>10.2f} mm³ ({nb:>5d}) {ve:>10.2f} mm³ ({ne:>5d}) "
          f"{ve/vb:>7.2f}×")

print("\n" + "=" * 76)
print("② tie 界面滑移：变形后「软骨底面节点 → 骨关节面」的距离")
print("=" * 76)
els, XYZ, sf = parse(W / "_E.feb")
NY = len(XYZ)
X0 = np.zeros((max(XYZ) + 1, 3))
for i, p in XYZ.items():
    X0[i] = p
X.to_hdf5(str(W / "_E.xplt"), str(W / "_rE2.h5"))
with h5py.File(W / "_rE2.h5", "r") as h:
    ns = sorted((int(q) for q in h["states"]), key=int)
    U0 = np.asarray(h[f"states/{ns[0]}/node_data/displacement/all"][()], float)
    U1 = np.asarray(h[f"states/{ns[-1]}/node_data/displacement/all"][()], float)


def pad(U):
    o = np.zeros((len(X0), 3))
    o[:len(U)] = U
    return o


U0, U1 = pad(U0), pad(U1)
for j in ("tibiotalar_tibia", "tibiotalar_talus"):
    pri = sf[f"tiesurf_{j}"]
    sec = sf[f"bonesurf_{j}"]
    pn = sorted({x for f in pri for x in f})
    sn = sorted({x for f in sec for x in f})
    B0 = X0[sn]
    B1 = X0[sn] + U1[sn]
    print(f"\n  【{j}】tie primary {len(pri)} 面/{len(pn)} 节点  "
          f"secondary {len(sec)} 面/{len(sn)} 节点")
    for tag, Xb, Uu in (("t=0", B0, None), ("t=1", B1, None)):
        P = X0[pn] + (U0[pn] if tag == "t=0" else U1[pn])
        d, _ = cKDTree(Xb).query(P)
        print(f"    {tag}: 软骨底面节点 → 骨面节点距离  min {d.min():.4f} "
              f"中位 {np.median(d):.4f} p95 {np.percentile(d,95):.4f} max {d.max():.4f} mm")
    # 界面相对滑移 = 变形后距离 − 变形前距离
    P0, P1 = X0[pn] + U0[pn], X0[pn] + U1[pn]
    d0, _ = cKDTree(B0).query(P0)
    d1, _ = cKDTree(B1).query(P1)
    slip = d1 - d0
    print(f"    **界面张开量(d1−d0)**：min {slip.min():+.4f} 中位 "
          f"{np.median(slip):+.4f} p95 {np.percentile(slip,95):+.4f} "
          f"max {slip.max():+.4f} mm")
    # 相对位移（滑移向量）
    tree = cKDTree(B0)
    _, idx = tree.query(P0)
    rel = (P1 - P0) - (B1[idx] - B0[idx])
    rn = np.linalg.norm(rel, axis=1)
    print(f"    **界面相对位移**：中位 {np.median(rn):.4f} p95 "
          f"{np.percentile(rn,95):.4f} max {rn.max():.4f} mm")
    print(f"    该软骨节点自身位移 |u| 中位 {np.median(np.linalg.norm(P1-P0,axis=1)):.4f} mm")
