"""末态真实接触斑：把 tibia 关节面节点投到 talus 关节面三角上，统计间隙。"""
import inspect
import re
import textwrap
from pathlib import Path

import h5py
import numpy as np
import pyfebio.xplt as X

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


def pt_tri(P, A, B, C):
    """P:(n,3) vs 三角 (A,B,C) 单三角 -> 距离 (n,)"""
    P = np.atleast_2d(np.asarray(P, float))
    AB, AC = B - A, C - A
    AP = P - A
    d1 = (AP @ AB)[:, None]
    d2 = (AP @ AC)[:, None]
    BP = P - B
    d3 = (BP @ AB)[:, None]
    d4 = (BP @ AC)[:, None]
    CP = P - C
    d5 = (CP @ AB)[:, None]
    d6 = (CP @ AC)[:, None]
    VC = np.cross(AB, AC)
    l2 = float(VC @ VC)
    if l2 < 1e-30:
        return np.linalg.norm(AP, axis=1)
    Wp = A + (d1 * AC - d2 * AB) / l2
    Up = B + (d3 * AC - d4 * AB) / l2
    Vp = C + (d5 * AC - d6 * AB) / l2
    d = np.full(len(P), np.inf)
    m = ((d1[:, 0] >= 0) & (d2[:, 0] >= 0) & (d1[:, 0] + d2[:, 0] <= l2))
    d[m] = np.linalg.norm(P[m] - Wp[m], axis=1)
    m = (d3[:, 0] <= 0) & (d4[:, 0] >= d3[:, 0])
    d[m] = np.linalg.norm(P[m] - B, axis=1)
    m = (d6[:, 0] >= 0) & (d5[:, 0] <= d6[:, 0])
    d[m] = np.linalg.norm(P[m] - C, axis=1)
    m = ~np.isfinite(d)
    if m.any():
        Pm = P[m]
        s = np.clip(((Pm - A) @ AB) / (AB @ AB), 0, 1)[:, None]
        t = np.clip(((Pm - A) @ AC) / (AC @ AC), 0, 1)[:, None]
        d[m] = np.minimum(
            np.linalg.norm(Pm - (A + s * AB), axis=1),
            np.linalg.norm(Pm - (A + t * AC), axis=1))
    return d


def surf_pts(sf, XYZ, U):
    return {i: XYZ[i] + U[i - 1] for i in sorted({x for f in sf for x in f})}


for tag, label in (("B", "基线 PENALTY=1 (F=419.56)"),
                   ("C", "基线 AUGLAG (F=315.06)")):
    txt = (W / f"_{tag}.feb").read_text(encoding="latin-1")
    XYZ = {}
    for m in re.finditer(
            r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)', txt):
        XYZ[int(m.group(1))] = np.array([float(m.group(2)), float(m.group(3)),
                                         float(m.group(4))])
    sf = {}
    for m in re.finditer(r'<Surface name="([^"]+)"[^>]*>(.*?)</Surface>', txt, re.S):
        sf[m.group(1)] = [tuple(int(y) for y in s.split(",")) for s in
                          re.findall(r"<[a-z0-9]+ id=[^>]*>\s*([\d,]+)\s*</", m.group(2))]
    X.to_hdf5(str(W / f"_{tag}.xplt"), str(W / f"_p{tag}.h5"))
    with h5py.File(W / f"_p{tag}.h5", "r") as h:
        ns = sorted((int(q) for q in h["states"]), key=int)
        u0 = np.asarray(h[f"states/{ns[0]}/node_data/displacement/all"][()], float)
        u1 = np.asarray(h[f"states/{ns[-1]}/node_data/displacement/all"][()], float)
    NX = max(XYZ) + 1
    U0 = np.zeros((NX, 3))
    U0[:len(u0)] = u0
    U1 = np.zeros((NX, 3))
    U1[:len(u1)] = u1

    plt_ = sf["surf_tibiotalar_tibia"]
    pls_ = sf["surf_tibiotalar_talus"]
    pn = sorted({x for f in plt_ for x in f})
    tri = np.array([np.array([XYZ[i] + U1[i - 1] for i in f]) for f in pls_])

    print("\n" + "=" * 78)
    print(f"{tag}: {label}")
    print("=" * 78)
    for khi, U in (("t=0", U0), ("t=1", U1)):
        P = np.array([XYZ[i] + U[i - 1] for i in pn])
        d = np.empty(len(P))
        for k, q in enumerate(P):
            d[k] = min(pt_tri(q[None, :], t[0], t[1], t[2])[0] for t in tri)
        print(f"  {khi}: tibia 关节面 {len(P)} 节点 → talus 关节面 距离 "
              f"min {d.min():.4f} 中位 {np.median(d):.4f} max {d.max():.4f} mm")
        for th in (0.01, 0.05, 0.1, 0.5, 1.0):
            k = int((d < th).sum())
            print(f"      < {th:<5.2f} mm : {k:>3d}/{len(d)} ({k/len(d)*100:4.1f}%)")
    # 面积估计
    A = 0.0
    for t in tri:
        A += 0.5 * np.linalg.norm(np.cross(t[1] - t[0], t[2] - t[0]))
    print(f"  talus 关节面总面积 = {A:.1f} mm²")
