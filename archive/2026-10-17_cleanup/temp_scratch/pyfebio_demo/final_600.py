"""跑到 600 N 档（prism），并算接触面积 + 软骨应力，完成双靶对照。

只用 Python 驱动，不经 .bat。
"""
import os
import re
import sys
import textwrap
import inspect
import subprocess
import h5py
import numpy as np
from pathlib import Path
from scipy.spatial import cKDTree

ROOT = Path(r"D:\Project\climbing_fall_analysis")
W = ROOT / "temp/thums/febio_lig"
PY = ROOT / ".venv/Scripts/python.exe"
MODEL = ROOT / "scripts/ankle_fe/ankle_lig_feb.py"
sys.path.insert(0, str(Path(__file__).parent))
from pt_tri import dist_point_to_tris

env = dict(os.environ)
env.update(dict(MATSET="WK", LIG_MODE="shell", DRIVER="disp",
                TOPFIX="lateral", TRUNC_MM="60", FIX_TALUS="all",
                CART_ELEM="prism", TIME_STEPS="30", DISP_MM="1.20"))
env["PYTHONIOENCODING"] = "utf-8"
r = subprocess.run([str(PY), "-u", str(MODEL)], cwd=str(ROOT), env=env,
                   capture_output=True, text=True, encoding="utf-8",
                   errors="replace", timeout=3600)
out = (r.stdout or "") + (r.stderr or "")
m = re.search(r"FEBio rc=(\d+)\s+用时 ([\d.]+)s", out)
print(f"FEBio rc={m.group(1)} 用时 {m.group(2)}s")
for l in out.split("\n"):
    if "tibiotalar" in l and "保留" in l:
        print("  ", l.strip())
    if "接触面" in l:
        print("  ", l.strip())

# ---------- 读末态 ----------
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
dst = W / "_fin.h5"
X.to_hdf5(str(src), str(dst))
AX = np.array([0.8153, 0.135, 0.5631])

with h5py.File(dst, "r") as h:
    st = sorted(int(k) for k in h["states"])
    last = max(st)
    D = np.asarray(h[f"states/{last}/node_data/displacement/all"][()], float)
    R = np.asarray(
        h[f"states/{last}/node_data/reaction forces/all"][()], float)
    D = D.reshape(-1, 3)
    F = float((R.reshape(-1, 3) @ AX).sum())
    print(f"\n状态 {len(st)} 末态 states/{last}")
    print(f"|u|max = {np.linalg.norm(D,axis=1).max():.4f} mm")
    print(f"支反力(沿加载轴) = {F:.1f} N")
    for k in h[f"states/{last}/element_data/stress"]:
        if "tibiotalar" in k:
            S = np.asarray(
                h[f"states/{last}/element_data/stress/{k}"][()], float)
            if S.shape[1] >= 6:
                sx, sy, sz, sxy, syz, szx = (S[:, 0], S[:, 1], S[:, 2],
                                             S[:, 3], S[:, 4], S[:, 5])
                vm = np.sqrt(0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2
                                    + (sz - sx) ** 2
                                    + 6 * (sxy ** 2 + syz ** 2 + szx ** 2)))
                print(f"  {k:30s} p99 {np.percentile(vm,99):6.3f}  "
                      f"max {vm.max():6.3f} MPa")

# ---------- 接触面积 ----------
txt = (W / "ankle_lig_WK.feb").read_text(encoding="utf-8")
nsec = re.search(r"<Nodes[^>]*>(.*?)</Nodes>", txt, re.S).group(1)
nid, co = [], []
for mm in re.finditer(r'<node id="(\d+)">([^<]+)</node>', nsec):
    v = [float(x) for x in re.split(r"[,\s]+", mm.group(2).strip()) if x]
    nid.append(int(mm.group(1)))
    co.append(v[:3])
nid = np.array(nid)
X0 = np.array(co)
lut = {n: k for k, n in enumerate(nid)}
SURF = {}
for nm in ("surf_tibiotalar_tibia", "surf_tibiotalar_talus"):
    blk = re.search(rf'<Surface name="{nm}">(.*?)</Surface>', txt, re.S)
    SURF[nm] = np.array([[int(x) for x in re.split(r"[,\s]+", e.group(1).strip())
                          if x]
                         for e in re.finditer(r"<tri3[^>]*>([^<]+)</tri3>",
                                              blk.group(1))
                         if len([int(x) for x in
                                 re.split(r"[,\s]+", e.group(1).strip()) if x]) == 3])

# 节点顺序必须与 xplt 一致
cp = None
with h5py.File(dst, "r") as h:
    nd = h["meshes/0/nodes"][()]
if getattr(nd, "dtype", None) is not None and nd.dtype.names:
    Xp = np.column_stack([nd["x"], nd["y"], nd["z"]]).astype(float)
else:
    Xp = np.asarray(nd, float).reshape(-1, 3)
print(f"xplt 节点 {Xp.shape}  .feb 节点 {len(nid)}")
if len(Xp) != len(nid):
    print("!! 节点数不一致，跳过接触面积")
    sys.exit(0)
Xd = Xp + D
AT = Xd[[[lut[n] for n in t] for t in SURF["surf_tibiotalar_tibia"]]]
AA = Xd[[[lut[n] for n in t] for t in SURF["surf_tibiotalar_talus"]]]


def area(T):
    return 0.5 * np.linalg.norm(np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]),
                                axis=1)


aT = area(AT)
tree = cKDTree(AA.mean(axis=1))
_, cand = tree.query(AT.mean(axis=1), k=min(40, len(AA)))
dmin = np.array([dist_point_to_tris(AT.mean(axis=1)[k],
                                    AA[cand[k]].reshape(-1, 3, 3)).min()
                 for k in range(len(AT))])
print("\n接触面积：")
for tol in (0.05, 0.10, 0.20, 0.50):
    mm_ = dmin < tol
    print(f"  间隙<{tol:5.2f} mm -> {mm_.sum():4d}/{len(dmin)} 面, "
          f"{aT[mm_].sum():7.1f} mm²")
print("\n吴恺靶子: 最大接触应力 3.74 MPa / 接触面积 335.5 mm²")
