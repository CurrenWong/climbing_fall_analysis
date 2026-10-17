"""算「胫骨下关节面接触面积」，对标吴恺 335.5 mm²（表2）。

读 .feb 拿两个 tibiotalar 接触面 + 节点，读 .xplt 拿末帧位移，
变形状下逐面求最近间隙，间隙 < tol 的面算接触。
"""
import re
import sys
import textwrap
import inspect
import h5py
import numpy as np
from pathlib import Path
from scipy.spatial import cKDTree

ROOT = Path(r"D:\Project\climbing_fall_analysis")
W = ROOT / "temp/thums/febio_lig"
tag = sys.argv[1] if len(sys.argv) > 1 else "WK"
FEB = W / f"ankle_lig_{tag}.feb"
if not FEB.exists():
    FEB = sorted(W.glob("ankle_lig_*.feb"),
                 key=lambda p: p.stat().st_mtime, reverse=True)[0]

# ---------- 1) 解析 .feb（逗号分隔） ----------
t = FEB.read_text(encoding="utf-8")
print(f"读 {FEB.name} ({FEB.stat().st_size / 1024:.0f} KB)")

nsec = re.search(r"<Nodes[^>]*>(.*?)</Nodes>", t, re.S).group(1)
nid, coords = [], []
for m in re.finditer(r'<node id="(\d+)">([^<]+)</node>', nsec):
    v = [float(x) for x in re.split(r"[,\s]+", m.group(2).strip()) if x]
    nid.append(int(m.group(1)))
    coords.append(v[:3])
nid = np.array(nid)
X0 = np.array(coords)
print(f"节点 {len(nid)}")

SURF = {}
for nm in ("surf_tibiotalar_tibia", "surf_tibiotalar_talus"):
    blk = re.search(rf'<Surface name="{nm}">(.*?)</Surface>', t, re.S)
    tris = []
    if blk:
        for m in re.finditer(r"<tri3[^>]*>([^<]+)</tri3>", blk.group(1)):
            ids = [int(x) for x in re.split(r"[,\s]+", m.group(1).strip()) if x]
            if len(ids) == 3:
                tris.append(ids)
    SURF[nm] = np.array(tris, dtype=int)
    print(f"  {nm}: {len(tris)} tri3")

# ---------- 2) 读 .xplt 位移 ----------
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

xplts = sorted(W.glob(f"ankle_lig_{tag}*.xplt"),
               key=lambda p: p.stat().st_mtime, reverse=True)
src = xplts[0]
dst = W / "_ca_result.h5"
print(f"读 {src.name} ({src.stat().st_size / 1024 / 1024:.1f} MB)")
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
if disp is None:
    print("!! 无位移；keys:", keys[:30])
    sys.exit(1)
D = np.asarray(disp, dtype=float)
if D.ndim == 3:
    D = D[-1]
D = D.reshape(-1, 3)
_c = np.asarray(coords_plt)
Xp = (np.column_stack([_c["x"], _c["y"], _c["z"]]).astype(float)
      if _c.dtype.names else _c.reshape(-1, 3).astype(float))
print(f"xplt 节点 {Xp.shape}  位移 {D.shape}")

# ---------- 3) 变形后表面 ----------
Xd = Xp + D
lut = {n: k for k, n in enumerate(nid)}
idxT = np.array([[lut[n] for n in tri] for tri in SURF["surf_tibiotalar_tibia"]])
idxA = np.array([[lut[n] for n in tri] for tri in SURF["surf_tibiotalar_talus"]])
AT = Xd[idxT]          # (F,3,3) 胫骨侧
AA = Xd[idxA]          # 距骨侧


def area(T):
    return 0.5 * np.linalg.norm(np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]),
                                axis=1)


aT = area(AT)
tree = cKDTree(AA.reshape(-1, 3))
d, _ = tree.query(AT.mean(axis=1))

print("\n间隙（胫骨软骨面质心 -> 距骨软骨面最近距离），mm")
for q in (0, 1, 5, 25, 50, 75, 100):
    print(f"  p{q:>3}: {np.percentile(d, q):>8.3f}")
print()
for tol in (0.02, 0.05, 0.10, 0.20, 0.50, 1.00):
    m = d < tol
    print(f"  间隙<{tol:5.2f} mm -> {m.sum():4d}/{len(d)} 面, "
          f"面积 {aT[m].sum():7.1f} mm²")
print(f"\n胫骨侧软骨面总面积 = {aT.sum():.1f} mm²")
print("吴恺 335.5 mm²")
