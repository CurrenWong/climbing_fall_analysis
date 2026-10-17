"""接触面积（精确点到三角面距离版），对标吴恺 335.5 mm²。"""
import re
import sys
import textwrap
import inspect
import h5py
import numpy as np
from pathlib import Path
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).parent))
from pt_tri import dist_point_to_tris

ROOT = Path(r"D:\Project\climbing_fall_analysis")
W = ROOT / "temp/thums/febio_lig"
FEB = sorted(W.glob("ankle_lig_*.feb"),
             key=lambda p: p.stat().st_mtime, reverse=True)[0]
t = FEB.read_text(encoding="utf-8")
print(f"读 {FEB.name}")

nsec = re.search(r"<Nodes[^>]*>(.*?)</Nodes>", t, re.S).group(1)
nid, coords = [], []
for m in re.finditer(r'<node id="(\d+)">([^<]+)</node>', nsec):
    v = [float(x) for x in re.split(r"[,\s]+", m.group(2).strip()) if x]
    nid.append(int(m.group(1)))
    coords.append(v[:3])
nid = np.array(nid)
X0 = np.array(coords)

SURF = {}
for nm in ("surf_tibiotalar_tibia", "surf_tibiotalar_talus"):
    blk = re.search(rf'<Surface name="{nm}">(.*?)</Surface>', t, re.S)
    tris = []
    for m in re.finditer(r"<tri3[^>]*>([^<]+)</tri3>", blk.group(1)):
        ids = [int(x) for x in re.split(r"[,\s]+", m.group(1).strip()) if x]
        if len(ids) == 3:
            tris.append(ids)
    SURF[nm] = np.array(tris, dtype=int)
print(f"  tibiotalar: {len(SURF['surf_tibiotalar_tibia'])} / "
      f"{len(SURF['surf_tibiotalar_talus'])} tri3")

# ---- 位移 ----
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

src = sorted(W.glob("ankle_lig_*.xplt"),
             key=lambda p: p.stat().st_mtime, reverse=True)[0]
dst = W / "_ca.h5"
print(f"读 {src.name} ({src.stat().st_size / 1024 / 1024:.1f} MB)")
X.to_hdf5(str(src), str(dst))

cplt = disp = None
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
            cplt = h[k][()]
        elif "displacement" in kl:
            disp = h[k][()]
D = np.asarray(disp, float)
if D.ndim == 3:
    D = D[-1]
D = D.reshape(-1, 3)
_c = np.asarray(cplt)
Xp = (np.column_stack([_c["x"], _c["y"], _c["z"]]).astype(float)
      if _c.dtype.names else _c.reshape(-1, 3).astype(float))
print(f"末帧 |u|max = {np.linalg.norm(D, axis=1).max():.3f} mm")

Xd = Xp + D
lut = {n: k for k, n in enumerate(nid)}
AT = Xd[[[lut[n] for n in tri] for tri in SURF["surf_tibiotalar_tibia"]]]
AA = Xd[[[lut[n] for n in tri] for tri in SURF["surf_tibiotalar_talus"]]]


def area(T):
    return 0.5 * np.linalg.norm(np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]),
                                axis=1)


aT = area(AT)
# 候选粗筛（用三角形质心）
tree = cKDTree(AA.mean(axis=1))
centT = AT.mean(axis=1)
_, cand = tree.query(centT, k=min(40, len(AA)))

print("\n精确点到三角面间隙（胫骨软骨面质心 -> 距骨软骨面）")
dmin = np.empty(len(centT))
for k in range(len(centT)):
    dmin[k] = dist_point_to_tris(centT[k],
                                AA[cand[k]].reshape(-1, 3, 3)).min()
for q in (0, 1, 5, 25, 50, 75, 100):
    print(f"  p{q:>3}: {np.percentile(dmin, q):>8.3f} mm")

print()
for tol in (0.02, 0.05, 0.10, 0.20, 0.50):
    m = dmin < tol
    print(f"  间隙<{tol:5.2f} mm -> {m.sum():4d}/{len(dmin)} 面, "
          f"面积 {aT[m].sum():7.1f} mm²")
print(f"\n胫骨侧软骨面总面积 = {aT.sum():.1f} mm²   吴恺 335.5 mm²")
