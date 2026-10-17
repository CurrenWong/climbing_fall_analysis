"""帧结构检查：D 的 shape / 每帧 |u|max。"""
import inspect
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

src = sorted(W.glob("ankle_tibiotalar.xplt"), key=lambda p: -p.stat().st_mtime)[0]
print(f"xplt {src.name}  {src.stat().st_size/1024/1024:.1f} MB")
dst = W / "_tie.h5"
X.to_hdf5(str(src), str(dst))
with h5py.File(dst, "r") as h:
    print("\nHDF5 树：")
    def walk(g, pre=""):
        for k in g:
            if isinstance(g[k], h5py.Group):
                print(f"  {pre}/{k}/")
                walk(g[k], pre + "/" + k)
            else:
                print(f"  {pre}/{k:32s} {g[k].shape} {g[k].dtype}")
    walk(h)
    ds = {}
    def _w(g, pre=""):
        for k in g:
            if isinstance(g[k], h5py.Group):
                _w(g[k], pre + "/" + k)
            else:
                try:
                    ds[pre + "/" + k] = np.asarray(g[k][()], dtype=float)
                except (TypeError, ValueError):
                    pass
    _w(h)
    for k, v in ds.items():
        if "displacement" in k.lower():
            print(f"\n{'='*70}\n位移数据集 {k}  shape={v.shape}")
            D = v
            if D.ndim == 3:
                for i in range(D.shape[0]):
                    print(f"  帧{i:3d} |u|max = {np.linalg.norm(D[i].reshape(-1,3),axis=1).max():.6f}")
            else:
                print(f"  单帧 |u|max = {np.linalg.norm(D.reshape(-1,3),axis=1).max():.6f}")
        if "reaction" in k.lower():
            print(f"反应力数据集 {k} shape={v.shape}")
