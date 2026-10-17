"""用已验证的 pt-tri 实现：前向/反向位移下关节面真实间隙演化。"""
import inspect
import re
import shutil
import textwrap
from pathlib import Path

import h5py
import numpy as np
import pyfebio.xplt as X

import pttri

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio")
AX = np.array([0.8153, 0.135, 0.5631])
AX = AX / np.linalg.norm(AX)

for ext in (".feb", ".xplt"):
    s = W / f"ankle_tibiotalar{ext}"
    if s.exists():
        shutil.copy2(s, W / f"_R{ext}")

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

for tag, label in (("B", "前向 DISP=+0.1472"), ("R", "反向 DISP=−0.1472"),
                   ("E", "m2 前向 TIE=1 AUGLAG")):
    fp, xp = W / f"_{tag}.feb", W / f"_{tag}.xplt"
    if not (fp.exists() and xp.exists()):
        print(f"{tag}: 缺文件")
        continue
    txt = fp.read_text(encoding="latin-1")
    XYZ = {}
    for m in re.finditer(
            r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)', txt):
        XYZ[int(m.group(1))] = np.array([float(m.group(2)), float(m.group(3)),
                                         float(m.group(4))])
    sf = {}
    for m in re.finditer(r'<Surface name="([^"]+)"[^>]*>(.*?)</Surface>', txt, re.S):
        sf[m.group(1)] = [tuple(int(y) for y in s.split(",")) for s in
                          re.findall(r"<[a-z0-9]+ id=[^>]*>\s*([\d,]+)\s*</", m.group(2))]
    X.to_hdf5(str(xp), str(W / f"_g{tag}.h5"))
    with h5py.File(W / f"_g{tag}.h5", "r") as h:
        ns = sorted((int(q) for q in h["states"]), key=int)
        u1 = np.asarray(h[f"states/{ns[-1]}/node_data/displacement/all"][()], float)
        R = np.asarray(h[f"states/{ns[-1]}/node_data/reaction forces/all"][()], float)
    NX = max(XYZ) + 1
    U1 = np.zeros((NX, 3))
    U1[:len(u1)] = u1
    pn = sorted({x for f in sf["surf_tibiotalar_tibia"] for x in f})
    tris = np.array([np.array([XYZ[i] + U1[i - 1] for i in f])
                     for f in sf["surf_tibiotalar_talus"]])
    P = np.array([XYZ[i] + U1[i - 1] for i in pn])
    d = pttri.dist_pt_tris(P, tris)
    print("\n" + "=" * 78)
    print(f"{tag}: {label}   states={len(ns)}   F={(R@AX).sum():.2f} N   "
          f"|u|max={np.linalg.norm(U1,axis=1).max():.4f}")
    print("=" * 78)
    print(f"  末态 关节面间隙(tibia {len(P)} 节点 → talus {len(tris)} 三角):  "
          f"min {d.min():.4f}  p5 {np.percentile(d,5):.4f} "
          f"中位 {np.median(d):.4f}  max {d.max():.4f} mm")
    for th in (0.001, 0.01, 0.05, 0.1, 0.5, 1.0):
        k = int((d < th).sum())
        print(f"      < {th:<6.3f} mm : {k:>3d}/{len(d)} ({k/len(d)*100:4.1f}%)")
