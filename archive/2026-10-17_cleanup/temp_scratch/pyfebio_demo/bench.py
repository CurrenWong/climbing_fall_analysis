"""基准脚本：跑多组配置，读末态「支反力 + 软骨应力 + 接触面数」。

不经过 .bat（cmd 的 GBK 编码会把含中文的批处理解析坏）。
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

ROOT = Path(r"D:\Project\climbing_fall_analysis")
W = ROOT / "temp/thums/febio_lig"
PY = ROOT / ".venv/Scripts/python.exe"
MODEL = ROOT / "scripts/ankle_fe/ankle_lig_feb.py"

BASE = dict(MATSET="WK", LIG_MODE="shell", DRIVER="disp", TOPFIX="lateral",
            TRUNC_MM="60", FIX_TALUS="all", TIME_STEPS="20")


def run(tag, **over):
    env = dict(os.environ)
    env.update(BASE)
    env.update({k: str(v) for k, v in over.items()})
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([str(PY), "-u", str(MODEL)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=3600)
    out = (r.stdout or "") + (r.stderr or "")
    rc = None
    m = re.search(r"FEBio rc=(\d+)\s+用时 ([\d.]+)s", out)
    if m:
        rc, secs = int(m.group(1)), float(m.group(2))
    else:
        secs = None
    carts = re.findall(r"(tibiotalar_\w+)\s+penta6\s+保留\s+(\d+)\s+剔除\s+(\d+)",
                       out)
    carts4 = re.findall(r"(tibiotalar_\w+)\s+tet4\s+保留\s+(\d+)\s+剔除\s+(\d+)",
                        out)
    surfs = re.findall(r"接触面 (tibiotalar_\w+): (\d+) 面", out)
    return dict(tag=tag, rc=rc, secs=secs, cars_keep=carts or carts4,
                surfs=surfs, out=out)


# ---- 读末态的 pyfebio 补丁（只做一次） ----
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

AX = np.array([0.8153, 0.135, 0.5631])


def vm(S):
    sx, sy, sz, sxy, syz, szx = (S[:, 0], S[:, 1], S[:, 2],
                                 S[:, 3], S[:, 4], S[:, 5])
    return np.sqrt(0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2
                          + 6 * (sxy ** 2 + syz ** 2 + szx ** 2)))


def read_last(tag):
    xplts = sorted(W.glob("ankle_lig_WK*.xplt"), key=lambda p: -p.stat().st_mtime)
    if not xplts:
        return {}
    dst = W / f"_b_{tag}.h5"
    X.to_hdf5(str(xplts[0]), str(dst))
    res = {}
    with h5py.File(dst, "r") as h:
        st = sorted(int(k) for k in h["states"])
        last = max(st)
        D = np.asarray(h[f"states/{last}/node_data/displacement/all"][()], float)
        R = np.asarray(
            h[f"states/{last}/node_data/reaction forces/all"][()], float)
        res["nstate"] = len(st)
        res["umax"] = float(np.linalg.norm(D.reshape(-1, 3), axis=1).max())
        res["F"] = float((R.reshape(-1, 3) @ AX).sum())
        if "element_data" in h[f"states/{last}"]:
            for k in h[f"states/{last}/element_data/stress"]:
                if "penta6" in k or "tet4" in k and "tibiotalar" in k:
                    S = np.asarray(
                        h[f"states/{last}/element_data/stress/{k}"][()], float)
                    if S.shape[1] >= 6 and "tibiotalar" in k:
                        v = vm(S)
                        res.setdefault("cart", {})[k] = (
                            float(np.percentile(v, 99)), float(v.max()))
                if "tibia_cort" in k:
                    S = np.asarray(
                        h[f"states/{last}/element_data/stress/{k}"][()], float)
                    if S.shape[1] >= 6:
                        res["tibia_cort_max"] = float(vm(S).max())
    return res


CONFIGS = [
    ("prism_dm0.5", dict(CART_ELEM="prism", DISP_MM="0.1472", DET_MIN="0.5")),
    ("prism_dm0.2", dict(CART_ELEM="prism", DISP_MM="0.1472", DET_MIN="0.2")),
    ("prism_dm0.1", dict(CART_ELEM="prism", DISP_MM="0.1472", DET_MIN="0.1")),
    ("prism_dm0.05", dict(CART_ELEM="prism", DISP_MM="0.1472", DET_MIN="0.05")),
    ("prism_dm0.02", dict(CART_ELEM="prism", DISP_MM="0.1472", DET_MIN="0.02")),
    ("tet4_dm0.5", dict(CART_ELEM="tet4", DISP_MM="0.1472", DET_MIN="0.5")),
]

print(f"{'配置':14s} {'rc':>3} {'秒':>6} {'F(N)':>8} {'|u|max':>8} "
      f"{'软骨p99':>9} {'软骨max':>9}")
print("-" * 78)
rows = []
for tag, ov in CONFIGS:
    r = run(tag, **ov)
    d = read_last(tag)
    c = d.get("cart", {})
    cp99 = max((v[0] for v in c.values()), default=float("nan"))
    cmax = max((v[1] for v in c.values()), default=float("nan"))
    print(f"{tag:14s} {str(r['rc']):>3} {str(r['secs']):>6} "
          f"{d.get('F', float('nan')):8.1f} {d.get('umax', float('nan')):8.3f} "
          f"{cp99:9.3f} {cmax:9.3f}")
    print(f"   {r['cars_keep']}  接触面 {r['surfs']}")
    rows.append((tag, r, d))
print()
print("对照：吴恺 最大接触应力 3.74 MPa、接触面积 335.5 mm²")
