"""读 FEBio 结果：末帧 |u|max、支反力沿加载轴分量、软骨 von Mises p99。"""
import textwrap
import inspect
import h5py
import numpy as np
from pathlib import Path

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio_lig")
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
import datetime
print(f"xplt: {src.name}  mtime "
      f"{datetime.datetime.fromtimestamp(src.stat().st_mtime):%H:%M:%S}  "
      f"{src.stat().st_size/1024/1024:.1f} MB")
dst = W / "_rr.h5"
X.to_hdf5(str(src), str(dst))

with h5py.File(dst, "r") as h:
    groups = {}

    def walk(g, pre=""):
        for k in g:
            if isinstance(g[k], h5py.Group):
                walk(g[k], pre + "/" + k)
            else:
                groups[pre + "/" + k] = g[k]
    walk(h)
    print("\n数据集：")
    for k, v in groups.items():
        print(f"  {k:60s} {v.shape}")

    def pick(tok):
        for k in groups:
            if tok in k.lower():
                return np.asarray(groups[k][()], dtype=float)
        return None

    umag = pick("displacement")
    reac = pick("reaction")
    # 用节点坐标定加载轴
    ax = np.array([0.8153, 0.135, 0.5631])

nstep = 1
print()
if umag is not None:
    A = umag.reshape(-1, 3) if umag.ndim != 3 else umag
    print(f"位移 shape {umag.shape}")
    last = umag[-1] if umag.ndim == 3 else umag
    print(f"末帧 |u|max = {np.linalg.norm(last.reshape(-1,3),axis=1).max():.4f} mm")
if reac is not None:
    print(f"支反力 shape {reac.shape}")
    last = reac[-1] if reac.ndim == 3 else reac
    last = last.reshape(-1, 3)
    s = (last @ ax).sum()
    print(f"支反力沿加载轴分量之和 = {s:.3f} N")
    print(f"  各分量 |x|={np.abs(last[:,0]).sum():.3f} "
          f"|y|={np.abs(last[:,1]).sum():.3f} |z|={np.abs(last[:,2]).sum():.3f}")
