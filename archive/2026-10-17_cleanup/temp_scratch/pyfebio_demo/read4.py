"""读 _B/_C/_D/_E 四个结果：F / |u|max / 轴向分量 / 高位移归属。"""
import inspect
import re
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

AX = np.array([0.8153, 0.135, 0.5631])
AX = AX / np.linalg.norm(AX)

CB = {"B": ("基线 PENALTY", "baseline"), "C": ("基线 AUGLAG", "baseline"),
      "D": ("m2 PENALTY", "m2"), "E": ("m2 AUGLAG", "m2")}
rows = []
for k, (label, mesh) in CB.items():
    xp = W / f"_{k}.xplt"
    fp = W / f"_{k}.feb"
    if not xp.exists():
        print(f"{k} {label}: 无 xplt（未收敛/无输出）")
        continue
    txt = fp.read_text(encoding="latin-1")
    parts = {}
    for m in re.finditer(
            r'<Elements\b[^>]*type="(tet4|hex8|penta6)"[^>]*name="([^"]+)"[^>]*>(.*?)</Elements>',
            txt, re.S):
        parts[m.group(2)] = {int(y) for s in
                             re.findall(r"<elem\b[^>]*>\s*([\d,]+)\s*</elem>", m.group(3))
                             for y in s.split(",")}
    X.to_hdf5(str(xp), str(W / f"_r{k}.h5"))
    with h5py.File(W / f"_r{k}.h5", "r") as h:
        ns = sorted((int(q) for q in h["states"]), key=int)
        fm = ns[-1]
        U = np.asarray(h[f"states/{fm}/node_data/displacement/all"][()], float)
        R = np.asarray(h[f"states/{fm}/node_data/reaction forces/all"][()], float)
    mag = np.linalg.norm(U, axis=1)
    uax = np.abs(U @ AX)
    reac = R @ AX
    print(f"\n{'='*76}\n{k}: {label}   states={len(ns)}  末帧={fm}\n{'='*76}")
    print(f"  F(沿轴支反力) = {reac.sum():.2f} N     "
          f"（+支反力之和 {reac[reac>0].sum():.1f} / −支反力之和 {reac[reac<0].sum():.1f}）")
    print(f"  |u|max = {mag.max():.4f} mm      |u·axis|max = {uax.max():.4f} mm")
    print(f"  胫骨侧 max = {max(mag[i-1] for i in parts.get('tibia_cort__hex8', [])):.4f}  "
          f"距骨侧 max = {max(mag[i-1] for i in parts.get('talus_cort__hex8', [])):.4f}")
    for dom in ("tibiotalar_tibia__tet4", "tibiotalar_talus__tet4"):
        if dom in parts:
            v = np.array([mag[i - 1] for i in sorted(parts[dom])])
            print(f"    {dom:26s} n={len(v):>5d} max={v.max():.4f} "
                  f"中位={np.median(v):.4f}")
    rows.append((k, label, mesh, reac.sum(), mag.max(), uax.max()))

print("\n\n" + "=" * 76)
print("汇总")
print("=" * 76)
print(f"  {'':3} {'配置':16s} {'网格':9s} {'F(N)':>9} {'|u|max':>9} {'|u·ax|max':>11}")
for r in rows:
    print(f"  {r[0]:>3} {r[1]:16s} {r[2]:9s} {r[3]:>9.2f} {r[4]:>9.4f} {r[5]:>11.4f}")
