"""末帧 |u|max 落在哪个部件 + 增量放大比（决定性判据）。"""
import inspect
import re
import textwrap
from pathlib import Path

import h5py
import numpy as np
import pyfebio.xplt as X

ROOT = Path(r"D:\Project\climbing_fall_analysis")
W = ROOT / "temp/thums/febio"
FEB = W / "ankle_tibiotalar.feb"

# ── 部件归属（从 .feb 的 Elements 块）──
txt = FEB.read_text(encoding="latin-1")
part_of = {}
for m in re.finditer(r'<Elements\b[^>]*type="(tet4|hex8|penta6)"[^>]*name="([^"]+)"[^>]*>(.*?)</Elements>',
                     txt, re.S):
    nm = m.group(2)
    for s in re.findall(r"<elem\b[^>]*>\s*([\d,]+)\s*</elem>", m.group(3)):
        for x in s.split(","):
            part_of[int(x)] = nm
# 软骨 offset 节点（tiesurf 里的）
carts = {x for n in ("tiesurf_tibiotalar_tibia", "tiesurf_tibiotalar_talus")
         for b in [re.search(rf'<Surface name="{n}"[^>]*>(.*?)</Surface>', txt, re.S).group(1)]
         for s in re.findall(r"<(?:tri3)\b[^>]*>\s*([\d,]+)\s*</", b)
         for x in [int(y) for y in s.split(",")]}
for x in carts:
    part_of.setdefault(x, "??cart??")

# ── xplt ──
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

src = sorted(W.glob("ankle_tibiotalar.xplt"), key=lambda p: -p.stat().st_mtime)[0]
X.to_hdf5(str(src), str(W / "_tie.h5"))
with h5py.File(W / "_tie.h5", "r") as h:
    ns = sorted((int(k) for k in h["states"]), key=int)
    U = {s: np.asarray(h[f"states/{s}/node_data/displacement/all"][()], dtype=float)
         for s in ns}
    Rn = {s: np.asarray(h[f"states/{s}/node_data/reaction forces/all"][()], dtype=float)
          for s in ns}
ax = np.array([0.8153, 0.135, 0.5631])
ax = ax / np.linalg.norm(ax)
f = ns[-1]

print("=" * 74)
print("① 增量放大比（施加位移 vs |u|max）  ⇒ 是否刚体模态")
print("=" * 74)
print(f"  {'帧':>4} {'t':>6} {'施加(mm)':>10} {'|u|max(mm)':>11} {'放大×':>8}")
for s in ns:
    t = s / (len(ns) - 1)
    given = t * 0.1472
    um = np.linalg.norm(U[s], axis=1).max()
    print(f"  {s:>4d} {t:>6.2f} {given:>10.5f} {um:>11.4f} "
          f"{um/given if given else float('nan'):>8.1f}")

mag = np.linalg.norm(U[f], axis=1)
print("\n" + "=" * 74)
print(f"② 末帧 |u|max = {mag.max():.4f} mm —— 落在哪个部件")
print("=" * 74)
order = np.argsort(-mag)[:15]
nid = np.asarray([int(x) for x in h["meshes/0/nodes"]["id"][()]]) if False else None
print(f"  {'idx':>7} {'|u|':>9}  u=(x,y,z)                     部件")
for i in order:
    print(f"  {i+1:>7d} {mag[i]:>9.4f}  ({U[f][i,0]:+.4f},{U[f][i,1]:+.4f},"
          f"{U[f][i,2]:+.4f})  {part_of.get(i+1,'?')}")

print("\n  各部位 |u| 统计（末帧）：")
bp = {}
for i, m in enumerate(mag):
    bp.setdefault(part_of.get(i + 1, "?"), []).append(m)
for nm, v in sorted(bp.items(), key=lambda kv: -max(kv[1])):
    v = np.array(v)
    print(f"    {nm:26s} n={len(v):>5d} max={v.max():>8.4f} "
          f"中位={np.median(v):>8.4f} p99={np.percentile(v,99):>8.4f}")

print("\n" + "=" * 74)
print("③ 支反力")
print("=" * 74)
R = Rn[f]
s_ax = (R @ ax).sum()
print(f"  末帧沿加载轴支反力 = {s_ax:.3f} N")
print(f"  分量和 |x|={np.abs(R[:,0]).sum():.2f} |y|={np.abs(R[:,1]).sum():.2f} "
      f"|z|={np.abs(R[:,2]).sum():.2f}")
# 按固定节点集合统计
print(f"  非零反力节点数 = {int((np.linalg.norm(R,axis=1)>1e-6).sum())}")
