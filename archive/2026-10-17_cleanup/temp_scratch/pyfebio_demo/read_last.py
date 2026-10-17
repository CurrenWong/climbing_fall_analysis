"""读末态：|u|max、支反力（沿加载轴）、软骨 von Mises p99。"""
import re
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
dst = W / "_r2.h5"
X.to_hdf5(str(src), str(dst))

import scipy.io
r = np.load(W / "lig_read_WK.npz", allow_pickle=True)
axis = np.asarray(r["axis"], dtype=float)

with h5py.File(dst, "r") as h:
    st = sorted(int(k) for k in h["states"])
    last = max(st)
    print(f"状态 {len(st)} 个，取末态 states/{last}")
    D = h[f"states/{last}/node_data/displacement/all"][()]
    R = h[f"states/{last}/node_data/reaction forces/all"][()]
    stg = {}
    if "element_data" in h[f"states/{last}"]:
        for k in h[f"states/{last}/element_data/stress"]:
            stg[k] = h[f"states/{last}/element_data/stress/{k}"][()]

D = np.asarray(D, float).reshape(-1, 3)
R = np.asarray(R, float).reshape(-1, 3)
print(f"末帧 |u|max = {np.linalg.norm(D, axis=1).max():.4f} mm")
print(f"支反力沿加载轴 = {(R @ axis).sum():.3f} N")

# 加载轴上的位移（近似轴向压缩）
print(f"R 各分量和: x {R[:,0].sum():.2f}  y {R[:,1].sum():.2f}  "
      f"z {R[:,2].sum():.2f}")

print("\n软骨 von Mises p99 / max (MPa)：")
for k, v in stg.items():
    if "cart" in k or "tibiotalar" in k:
        S = np.asarray(v, float)
        if S.shape[1] < 6:
            continue
        sx, sy, sz, sxy, syz, szx = (S[:, 0], S[:, 1], S[:, 2],
                                     S[:, 3], S[:, 4], S[:, 5])
        vm = np.sqrt(0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2
                            + (sz - sx) ** 2
                            + 6 * (sxy ** 2 + syz ** 2 + szx ** 2)))
        print(f"  {k:32s} n={len(vm):5d}  p99 {np.percentile(vm,99):7.3f}  "
              f"max {vm.max():7.3f}")
print("\n其它域的最大 von Mises（前 5）：")
allv = []
for k, v in stg.items():
    S = np.asarray(v, float)
    if S.shape[1] < 6:
        continue
    sx, sy, sz, sxy, syz, szx = (S[:, 0], S[:, 1], S[:, 2],
                                 S[:, 3], S[:, 4], S[:, 5])
    vm = np.sqrt(0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2
                        + 6 * (sxy ** 2 + syz ** 2 + szx ** 2)))
    allv.append((vm.max(), k))
for m, k in sorted(allv, reverse=True)[:6]:
    print(f"  {k:34s} max {m:9.3f}")
