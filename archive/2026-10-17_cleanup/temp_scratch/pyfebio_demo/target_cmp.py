"""双靶稳健性检验：让 tet4 与 penta6 各自跑到 ~600 N，比较【应力】与【接触面积】。

要点：刚度差 2.2 倍只影响"需要多大位移"；**靶子是否也差 2 倍**才决定
之前的验证结论能不能用。
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

AX = np.array([0.8153, 0.135, 0.5631])

# ---------- pyfebio 补丁 ----------
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


def vm_of(S):
    sx, sy, sz, sxy, syz, szx = (S[:, 0], S[:, 1], S[:, 2],
                                 S[:, 3], S[:, 4], S[:, 5])
    return np.sqrt(0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2
                          + 6 * (sxy ** 2 + syz ** 2 + szx ** 2)))


def run(elem, disp, steps=30):
    env = dict(os.environ)
    env.update(dict(MATSET="WK", LIG_MODE="shell", DRIVER="disp",
                    TOPFIX="lateral", TRUNC_MM="60", FIX_TALUS="all",
                    TIME_STEPS=str(steps), DISP_MM=str(disp),
                    CART_ELEM=elem))
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([str(PY), "-u", str(MODEL)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=3600)
    out = (r.stdout or "") + (r.stderr or "")
    m = re.search(r"FEBio rc=(\d+)", out)
    return (int(m.group(1)) if m else None), out


def measure(elem):
    src = sorted(W.glob("ankle_lig_WK*.xplt"),
                 key=lambda p: -p.stat().st_mtime)[0]
    dst = W / "_tb.h5"
    X.to_hdf5(str(src), str(dst))
    res = {}
    with h5py.File(dst, "r") as h:
        st = sorted(int(k) for k in h["states"])
        last = max(st)
        D = np.asarray(h[f"states/{last}/node_data/displacement/all"][()],
                       float).reshape(-1, 3)
        R = np.asarray(
            h[f"states/{last}/node_data/reaction forces/all"][()],
            float).reshape(-1, 3)
        res["F"] = float((R @ AX).sum())
        res["umax"] = float(np.linalg.norm(D, axis=1).max())
        best = 0.0
        for k in h[f"states/{last}/element_data/stress"]:
            if "tibiotalar" in k:
                S = np.asarray(
                    h[f"states/{last}/element_data/stress/{k}"][()], float)
                if S.shape[1] >= 6:
                    best = max(best, float(vm_of(S).max()))
        res["vm_max"] = best
        nd = h["meshes/0/nodes"][()]
    Xp = (np.column_stack([nd["x"], nd["y"], nd["z"]]).astype(float)
          if nd.dtype.names else np.asarray(nd, float).reshape(-1, 3))
    Xd = Xp + D

    txt = (W / "ankle_lig_WK.feb").read_text(encoding="utf-8")
    nsec = re.search(r"<Nodes[^>]*>(.*?)</Nodes>", txt, re.S).group(1)
    nid = [int(m.group(1))
           for m in re.finditer(r'<node id="(\d+)">', nsec)]
    lut = {n: k for k, n in enumerate(nid)}
    SURF = {}
    for nm in ("surf_tibiotalar_tibia", "surf_tibiotalar_talus"):
        blk = re.search(rf'<Surface name="{nm}">(.*?)</Surface>', txt, re.S)
        tris = []
        for e in re.finditer(r"<tri3[^>]*>([^<]+)</tri3>", blk.group(1)):
            ids = [int(x) for x in re.split(r"[,\s]+", e.group(1).strip()) if x]
            if len(ids) == 3:
                tris.append(ids)
        SURF[nm] = np.array(tris)
    AT = Xd[[[lut[n] for n in t] for t in SURF["surf_tibiotalar_tibia"]]]
    AA = Xd[[[lut[n] for n in t] for t in SURF["surf_tibiotalar_talus"]]]
    aT = 0.5 * np.linalg.norm(np.cross(AT[:, 1] - AT[:, 0],
                                       AT[:, 2] - AT[:, 0]), axis=1)
    tree = cKDTree(AA.mean(axis=1))
    _, cand = tree.query(AT.mean(axis=1), k=min(40, len(AA)))
    dm = np.array([dist_point_to_tris(AT.mean(axis=1)[k],
                                      AA[cand[k]].reshape(-1, 3, 3)).min()
                   for k in range(len(AT))])
    res["area020"] = float(aT[dm < 0.20].sum())
    res["area050"] = float(aT[dm < 0.50].sum())
    res["nface"] = (len(SURF["surf_tibiotalar_tibia"]),
                    len(SURF["surf_tibiotalar_talus"]))
    return res


print("=== 先给 tet4 找 ~600 N 的位移（它刚度约 2.2 倍）===")
for d in (0.30, 0.40):
    rc, out = run("tet4", d)
    if rc == 0:
        m = measure("tet4")
        print(f"  tet4 DISP={d:<5} -> F={m['F']:7.1f} N  |u|max={m['umax']:.3f}")
print("\n=== 双靶对照 ===")
print(f"{'单元':8s} {'DISP':>6} {'F(N)':>8} {'|u|max':>8} {'软骨vm max':>11} "
      f"{'面积@0.20':>11} {'面积@0.50':>11} {'接触面':>11}")
print("-" * 92)
for elem, d in (("tet4", 0.40), ("prism", 1.20)):
    rc, out = run(elem, d)
    if rc != 0:
        print(f"{elem:8s} {d:6.2f} rc={rc} 失败")
        continue
    m = measure(elem)
    print(f"{elem:8s} {d:6.2f} {m['F']:8.1f} {m['umax']:8.3f} "
          f"{m['vm_max']:11.3f} {m['area020']:11.1f} {m['area050']:11.1f} "
          f"{str(m['nface']):>11}")
print("\n吴恺靶子: 最大接触应力 3.74 MPa / 接触面积 335.5 mm²")
