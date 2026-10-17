"""加密网格 vs 原始网格：受控对比（同参数，仅网格分辨率不同，且可切单元类型）。"""
import os
import re
import textwrap
import inspect
import subprocess
import h5py
import numpy as np
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
MODEL = ROOT / "scripts/ankle_fe/ankle_lig_feb.py"
PY = ROOT / ".venv/Scripts/python.exe"
W = ROOT / "temp/thums/febio_lig"
AX = np.array([0.8153, 0.135, 0.5631])

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

BASE = dict(MATSET="WK", LIG_MODE="shell", DRIVER="disp", TOPFIX="lateral",
            TRUNC_MM="60", FIX_TALUS="all", TIME_STEPS="20", DISP_MM="0.1472",
            LIG_ONLY="1")


def run(label, elem, npz=None, cart=None):
    env = dict(os.environ); env.update(BASE)
    env["CART_ELEM"] = elem; env["PYTHONIOENCODING"] = "utf-8"
    # ★ 无条件清除/设置：否则继承自 shell 的 THUMS_NPZ 会串到"原始"组
    #   （已踩：四个 run 全用了 r1，连原始组也失败）
    env.pop("THUMS_NPZ", None)
    env.pop("CART_DIR", None)
    env.pop("CART_OUT", None)
    if npz:
        env["THUMS_NPZ"] = npz
    if cart:
        env["CART_DIR"] = cart
    r = subprocess.run([str(PY), "-u", str(MODEL)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=7200)
    o = (r.stdout or "") + (r.stderr or "")
    (ROOT / f"temp/pyfebio_demo/_res_{label}_{elem}.txt").write_text(o, encoding="utf-8")
    m = re.search(r"FEBio rc=(\d+)", o)
    rc = int(m.group(1)) if m else None
    t = re.search(r"用时 ([\d.]+)s", o)
    carts = re.findall(r"(tibiotalar_\w+)\s+(hex8|penta6|penta15|tet4)\s+保留\s+(\d+)\s+剔除\s+(\d+)", o)
    surfs = re.findall(r"接触面 (\w+): (\d+) 面", o)
    err = "负Jacobian" if "Negative jacobian" in o else (
        "不收敛" if "failed to converge" in o else "")
    F = um = None
    if rc == 0:
        src = sorted(W.glob("ankle_lig_WK*.xplt"), key=lambda p: -p.stat().st_mtime)[0]
        dst = W / "_cmp.h5"
        X.to_hdf5(str(src), str(dst))
        with h5py.File(dst, "r") as h:
            last = max(int(k) for k in h["states"])
            R = np.asarray(h[f"states/{last}/node_data/reaction forces/all"][()], float)
            Ud = np.asarray(h[f"states/{last}/node_data/displacement/all"][()], float)
            F = float((R.reshape(-1, 3) @ AX).sum())
            um = float(np.linalg.norm(Ud.reshape(-1, 3), axis=1).max())
    return rc, (t.group(1) if t else "?"), F, um, carts, surfs, err


R1N = str(ROOT / "temp/thums/thums_lowerlimb_r1.npz")
R1C = str(ROOT / "temp/thums/r1")

print("加密网格 vs 原始网格（DISP=0.1472 / TRUNC=60 / FIX_TALUS=all）")
print("=" * 92)
for tag, npz, cart in (("原始", None, None), ("加密R1", R1N, R1C)):
    for elem in ("tet4", "prism"):
        rc, t, F, um, carts, surfs, err = run(tag, elem, npz, cart)
        sF = f"{F:9.1f}" if F is not None else f"{'-':>9s}"
        print(f"{tag:>8s} {elem:>6s} rc={rc} {t:>7s}s  F={sF} N  "
              f"|u|max={um if um else '-'}  {err}")
        print(f"           软骨 {[(c[0][-5:], c[1], c[2], c[3]) for c in carts]}")
        print(f"           接触面 {surfs}")
    print()
