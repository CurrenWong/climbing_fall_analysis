"""收敛研究（AUGLAG）：CP_M = 1,2,3(+4) 扫，提取支反力/软骨应力/接触面积。

★ LAUGON=AUGLAG 是必须的：m≥2 时 PENALTY 在第一步就报 Negative jacobian
  （失败点恒为 1/(11·TIME_STEPS)，即近乎零载荷 ⇒ 不是几何问题，是接触求解器）。
"""
import inspect
import os
import re
import shutil
import subprocess
import textwrap
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
THUMS = ROOT / "temp/thums"
WORK = THUMS / "febio"
PY = str(ROOT / ".venv/Scripts/python.exe")

# ── pyfebio 兼容补丁（同 read_react.py）──
import pyfebio.xplt as X  # noqa: E402

_src = textwrap.dedent(inspect.getsource(X.parse_mesh))
for _o, _n in (('domain["name"]', 'domain.get("name", str(domain["id"][0]))'),
               ('surface["name"]', 'surface.get("name", str(surface["id"][0]))'),
               ("{domain['name']}", "{domain.get('name', str(domain['id'][0]))}"),
               ("{surface['name']}", "{surface.get('name', str(surface['id'][0]))}"),
               ('if nodeset["name"] == "":', 'if nodeset.get("name", "") == "":'),
               ("{nodeset['name']}",
                "{nodeset.get('name', str(nodeset['id'][0]))}"),
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

import h5py  # noqa: E402


def post(tag):
    import datetime
    src = WORK / "ankle_tibiotalar.xplt"
    dst = WORK / "_conv.h5"
    X.to_hdf5(str(src), str(dst))
    out = {"tag": tag}
    with h5py.File(dst, "r") as h:
        g = {}

        def walk(gr, pre=""):
            for k in gr:
                if isinstance(gr[k], h5py.Group):
                    walk(gr[k], pre + "/" + k)
                else:
                    g[pre + "/" + k] = gr[k]
        walk(h)

        def pick(tok, exclude=()):
            for k in sorted(g):
                kl = k.lower()
                if tok in kl and not any(e in kl for e in exclude):
                    return np.asarray(g[k][()], float), k
            return None, None

        um, _ = pick("displacement")
        re, _ = pick("reaction")
        if um is not None:
            last = um[-1] if um.ndim == 3 else um
            out["umax"] = float(np.linalg.norm(last.reshape(-1, 3), axis=1).max())
        if re is not None:
            last = re[-1] if re.ndim == 3 else re
            v = last.reshape(-1, 3).sum(axis=0)
            out["F_norm"] = float(np.linalg.norm(v))
            out["F_xyz"] = v
        # 软骨 von Mises p99：找包含 cartilage 的单元数据集
        for k in sorted(g):
            kl = k.lower()
            if "stress" in kl and "solid" in kl:
                A = np.asarray(g[k][()], float)
                out["stress_key"] = k
                out["stress_shape"] = str(A.shape)
                break
    return out


env0 = dict(os.environ)
env0["PATH"] = r"D:\Program\FEBioStudio\bin;" + env0["PATH"]
env0.update(CP_DOMAIN="bary", CP_DIAG="02", TIE="1", DRIVER="disp",
            DISP_MM="0.1472", TIME_STEPS="20", LAUGON="AUGLAG")

print(f"{'m':>3s} {'底三角':>8s} {'接触面':>12s} {'rc':>3s} {'耗时':>9s} "
      f"{'|u|max':>9s} {'|F|(N)':>9s}")
rows = []
for m in (1, 2, 3):
    out = f"cartilage_mesh_patch_m{m}.npz"
    e = dict(env0, CP_M=str(m), CP_OUT=out)
    r = subprocess.run([PY, "-u", str(ROOT / "scripts/ankle_fe/cart_patch.py")],
                       cwd=ROOT, env=e, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=3600)
    ncell = [int(x) for x in re.findall(r"面\s+(\d+)\(", r.stdout)]

    e2 = dict(env0, CART_MESH=out)
    r2 = subprocess.run([PY, "-u", str(ROOT / "scripts/ankle_fe/build_feb.py")],
                        cwd=ROOT, env=e2, capture_output=True, text=True,
                        encoding="utf-8", errors="replace", timeout=7200)
    o = r2.stdout + r2.stderr
    cf = re.findall(r"外层接触面 (\d+)", o)
    rc = re.search(r"FEBio rc=(\d+)", o)
    tm = re.search(r"用时 ([\d.]+)s", o)
    if rc and rc.group(1) == "0":
        p = post(f"m{m}")
        rows.append((m, sum(ncell), cf, rc.group(1), tm.group(1), p))
        print(f"{m:>3d} {sum(ncell):>8d} {str(cf):>12s} {rc.group(1):>3s} "
              f"{tm.group(1)+'s':>9s} {p.get('umax', 0):>9.4f} "
              f"{p.get('F_norm', 0):>9.3f}")
    else:
        bad = [l.strip() for l in o.split("\n")
               if "Negative" in l or "failed to converge" in l]
        print(f"{m:>3d} {sum(ncell):>8d} {str(cf):>12s} {rc.group(1) if rc else '?':>3s} "
              f"{(tm.group(1)+'s') if tm else '?':>9s}   FAIL {bad[:1]}")

print("\n--- 后处理探查 ---")
if rows:
    p = rows[-1][5]
    print("stress 数据集:", p.get("stress_key"), p.get("stress_shape"))
    print("F 分量:", p.get("F_xyz"))
    print("|F| =", p.get("F_norm"), " |u|max =", p.get("umax"))
