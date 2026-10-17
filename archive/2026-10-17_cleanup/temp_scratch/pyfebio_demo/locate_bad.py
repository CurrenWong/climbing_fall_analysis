"""定位 fibula_cort 里翻反的那 2 个 hex 单元（相对腓骨/截断面的位置）。"""
import re
import textwrap
import inspect
import h5py
import numpy as np
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
W = ROOT / "temp/thums/febio_lig"
t = (W / "ankle_lig_WK.feb").read_text(encoding="utf-8")

nsec = re.search(r"<Nodes[^>]*>(.*?)</Nodes>", t, re.S).group(1)
nid, coords = [], []
for m in re.finditer(r'<node id="(\d+)">([^<]+)</node>', nsec):
    v = [float(x) for x in re.split(r"[,\s]+", m.group(2).strip()) if x]
    nid.append(int(m.group(1)))
    coords.append(v[:3])
nid = np.array(nid)
X0 = np.array(coords)
lut = {n: k for k, n in enumerate(nid)}

DOM = {}
for m in re.finditer(r'<Elements type="(\w+)" name="([^"]+)">(.*?)</Elements>',
                     t, re.S):
    et, nm, body = m.group(1), m.group(2), m.group(3)
    els = [[int(x) for x in re.split(r"[,\s]+", e.group(1).strip()) if x]
           for e in re.finditer(r"<elem id=\"\d+\">([^<]+)</elem>", body)]
    DOM[nm] = (et, np.array(els, dtype=int))

# 位移
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
src = sorted(W.glob("ankle_lig_WK*.xplt"), key=lambda p: -p.stat().st_mtime)[0]
dst = W / "_dg2.h5"
X.to_hdf5(str(src), str(dst))
cp = disp = None
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
            cp = h[k][()]
        elif "displacement" in kl:
            disp = h[k][()]
D = np.asarray(disp, float)
if D.ndim == 3:
    D = D[-1]
D = D.reshape(-1, 3)
_c = np.asarray(cp)
Xp = (np.column_stack([_c["x"], _c["y"], _c["z"]]).astype(float)
      if _c.dtype.names else _c.reshape(-1, 3).astype(float))
Xd = Xp + D

et, els = DOM["fibula_cort__hex8"]
idx = np.array([[lut[n] for n in e] for e in els])
x = X0[idx]
y = Xd[idx]
pairs = [(0, 1), (2, 3), (4, 7)]
F = np.stack([y[:, b] - y[:, a] for a, b in pairs], 1)
G = np.stack([x[:, b] - x[:, a] for a, b in pairs], 1)
d = np.linalg.det(F) / np.linalg.det(G)

print("fibula_cort hex 单元")
print(f"  总数 {len(d)}，det(F) min {d.min():.4f}，"
      f"倒数第二 {np.sort(d)[1]:.4f}")
bad = np.argsort(d)[:4]
print(f"\n最差 4 个：")
for k in bad:
    cen = x[k].mean(axis=0)
    cen_d = y[k].mean(axis=0)
    print(f"  elem#{els[k][0]:<6d} det={d[k]:9.4f}  原始质心 "
          f"({cen[0]:8.2f},{cen[1]:8.2f},{cen[2]:8.2f})  "
          f"位移 {np.linalg.norm(cen_d - cen):.4f} mm")

allx = X0[np.array([[lut[n] for n in e] for e in els]).ravel()]
print(f"\n腓骨单元 y 范围 {allx[:, 1].min():.2f} ~ {allx[:, 1].max():.2f}")
for k in bad:
    print(f"  elem#{els[k][0]:<6d} 的 8 节点 y: "
          f"{np.round(x[k][:, 1], 2)}")
print(f"\n该单元 8 节点原始坐标：")
for k in bad[:2]:
    for j in range(8):
        print(f"   n{els[k][j]:<5d} {x[k][j]}  -> {y[k][j]}  "
              f"u {np.linalg.norm(y[k][j] - x[k][j]):.4f}")
    print()
