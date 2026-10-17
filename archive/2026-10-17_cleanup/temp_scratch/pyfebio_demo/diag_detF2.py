"""诊断：逐单元算 det(F)，找濒临/已翻转的单元（定位到域）。

.feb 结构：<Elements type="tet4" name="域名"> <elem id="1">a,b,c,d</elem> ... </Elements>
"""
import re
import sys
import textwrap
import inspect
import h5py
import numpy as np
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
W = ROOT / "temp/thums/febio_lig"
FEB = W / "ankle_lig_WK.feb"
t = FEB.read_text(encoding="utf-8")

nsec = re.search(r"<Nodes[^>]*>(.*?)</Nodes>", t, re.S).group(1)
nid, coords = [], []
for m in re.finditer(r'<node id="(\d+)">([^<]+)</node>', nsec):
    v = [float(x) for x in re.split(r"[,\s]+", m.group(2).strip()) if x]
    nid.append(int(m.group(1)))
    coords.append(v[:3])
nid = np.array(nid)
X0 = np.array(coords)
lut = {n: k for k, n in enumerate(nid)}
print(f"节点 {len(nid)}")

DOM = {}
for m in re.finditer(r'<Elements type="(\w+)" name="([^"]+)">(.*?)</Elements>',
                     t, re.S):
    et, nm, body = m.group(1), m.group(2), m.group(3)
    els = []
    for e in re.finditer(r"<elem id=\"\d+\">([^<]+)</elem>", body):
        ids = [int(x) for x in re.split(r"[,\s]+", e.group(1).strip()) if x]
        els.append(ids)
    DOM[nm] = (et, np.array(els, dtype=int))
    print(f"  {nm:34s} {et:6s} {len(els):6d}")

import pyfebio.xplt as X
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

xplts = sorted(W.glob("ankle_lig_WK*.xplt"), key=lambda p: -p.stat().st_mtime)
src = xplts[0]
dst = W / "_dg.h5"
print(f"\n读 {src.name} ({src.stat().st_size / 1024 / 1024:.1f} MB)")
X.to_hdf5(str(src), str(dst))
coords_plt = disp = None
with h5py.File(dst, "r") as h:
    keys = []

    def walk(g, pre=""):
        for k in g:
            if isinstance(g[k], h5py.Group):
                walk(g[k], pre + "/" + k)
            else:
                keys.append(pre + "/" + k)
    walk(h)
    for k in keys:
        kl = k.lower()
        if kl.endswith("/nodes"):
            coords_plt = h[k][()]
        elif "displacement" in kl:
            disp = h[k][()]
D = np.asarray(disp, float)
if D.ndim == 3:
    D = D[-1]
D = D.reshape(-1, 3)
_c = np.asarray(coords_plt)
Xp = (np.column_stack([_c["x"], _c["y"], _c["z"]]).astype(float)
      if _c.dtype.names else _c.reshape(-1, 3).astype(float))
Xd = Xp + D
print(f"末帧 |u|max = {np.linalg.norm(D, axis=1).max():.4f} mm")


def detF(idx, et):
    x = X0[idx]
    y = Xd[idx]
    if et == "tet4":
        F = np.stack([y[:, 1] - y[:, 0], y[:, 2] - y[:, 0],
                      y[:, 3] - y[:, 0]], 1)
        G = np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0],
                      x[:, 3] - x[:, 0]], 1)
    else:                                  # hex8
        pairs = [(0, 1), (2, 3), (4, 7)]
        F = np.stack([y[:, b] - y[:, a] for a, b in pairs], 1)
        G = np.stack([x[:, b] - x[:, a] for a, b in pairs], 1)
    return np.linalg.det(F) / np.linalg.det(G)


print("\n" + "=" * 88)
print(f"{'域':32s} {'类型':>6} {'单元':>6} {'det min':>11} {'p1':>10} "
      f"{'p5':>10} {'负个数':>7}")
print("-" * 88)
rows = []
for nm, (et, els) in DOM.items():
    if et not in ("tet4", "hex8"):
        continue
    idx = np.array([[lut[n] for n in ids] for ids in els])
    d = detF(idx, et)
    rows.append((d.min(), nm, et, len(d), float(np.percentile(d, 1)),
                 float(np.percentile(d, 5)), int((d <= 0).sum())))
for mn, nm, et, n, p1, p5, nneg in sorted(rows):
    print(f"{nm:32s} {et:>6} {n:6d} {mn:11.5f} {p1:10.5f} {p5:10.5f} "
          f"{nneg:7d}")
