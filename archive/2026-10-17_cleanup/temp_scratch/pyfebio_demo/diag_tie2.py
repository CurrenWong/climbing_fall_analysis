"""TIE 诊断 ②：tie 块 + surface 规模 + **找出 |u|max 落在哪个部件**（决定性）。"""
import re
from pathlib import Path

import numpy as np

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio")
FEB = W / "ankle_tibiotalar.feb"
txt = FEB.read_text(encoding="latin-1")

print("=" * 74)
print("① tie contact 块（原样）")
print("=" * 74)
for m in re.finditer(r'<contact\b[^>]*type="tied-elastic".*?</contact>', txt, re.S):
    print("  " + re.sub(r"\n\s*", "\n  ", m.group(0).strip()))

print("\n" + "=" * 74)
print("② SurfacePair + 各 Surface 规模")
print("=" * 74)
for m in re.finditer(r'<SurfacePair\b[^>]*name="([^"]+)"[^>]*>(.*?)</SurfacePair>', txt, re.S):
    a = re.search(r"<primary>([^<]+)", m.group(2))
    b = re.search(r"<secondary>([^<]+)", m.group(2))
    print(f"  {m.group(1):24s} primary={a.group(1) if a else '?':26s} "
          f"secondary={b.group(1) if b else '?'}")

surf = {}
for n, b in re.findall(r'<Surface name="([^"]+)"[^>]*>(.*?)</Surface>', txt, re.S):
    fs = re.findall(r"<(tri3|quad4|quad8)\b[^>]*>\s*([\d,]+)\s*</", b)
    surf[n] = ([tuple(int(y) for y in f.split(",")) for _, f in fs], [t for t, _ in fs])
for n in sorted(surf):
    fs, ts = surf[n]
    nn = len({x for f in fs for x in f})
    kinds = {}
    for t in ts:
        kinds[t] = kinds.get(t, 0) + 1
    print(f"  {n:26s} 面 {len(fs):>5d} {str(kinds):34s} 节点 {nn:>5d}")

# ── ③ 部件归属 ──
print("\n" + "=" * 74)
print("③ 部件 → 节点集")
print("=" * 74)
part_of = {}
for m in re.finditer(r'<Elements\b[^>]*type="(tet4|hex8|penta6)"[^>]*name="([^"]+)"[^>]*>(.*?)</Elements>',
                     txt, re.S):
    et, nm, body = m.group(1), m.group(2), m.group(3)
    ids = [tuple(int(y) for y in s.split(",")) for s in
           re.findall(r"<elem\b[^>]*>\s*([\d,]+)\s*</elem>", body)]
    n_el = len(ids)
    for e in ids:
        for x in e:
            part_of.setdefault(x, nm)
    print(f"  {nm:26s} {et:7s} 单元 {n_el:>7d}")

# ── ④ 读 xplt 末帧，找 |u|max 的节点 ──
print("\n" + "=" * 74)
print("④ |u|max 落在哪个部件（决定性）")
print("=" * 74)
import h5py
import textwrap
import inspect
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

src = sorted(W.glob("ankle_tibiotalar.xplt"), key=lambda p: -p.stat().st_mtime)[0]
print(f"  xplt: {src.name} {src.stat().st_size/1024/1024:.1f} MB")
dst = W / "_tie.h5"
X.to_hdf5(str(src), str(dst))
with h5py.File(dst, "r") as h:
    ds = {}

    def _walk(g, pre=""):
        for k in g:
            if isinstance(g[k], h5py.Group):
                _walk(g[k], pre + "/" + k)
            else:
                try:
                    ds[pre + "/" + k] = np.asarray(g[k][()], dtype=float)
                except (TypeError, ValueError):
                    pass
    _walk(h)
    D = next((v for k, v in ds.items() if "displacement" in k.lower()), None)
ax = np.array([0.8153, 0.135, 0.5631]) / np.linalg.norm([0.8153, 0.135, 0.5631])
assert D is not None, "无 displacement"
last = D[-1].reshape(-1, 3) if D.ndim == 3 else D.reshape(-1, 3)
mag = np.linalg.norm(last, axis=1)
nid = np.argsort(-mag)[:12] + 1          # 1-based 节点号
print(f"\n  末帧 |u|max = {mag.max():.4f} mm")
print("  位移最大的 12 个节点：")
for i in nid:
    print(f"    node {i:>6d}  |u|={mag[i-1]:.4f}  u=({last[i-1,0]:+.4f},"
          f"{last[i-1,1]:+.4f},{last[i-1,2]:+.4f})  {part_of.get(i, '?')}")
print("\n  各部件 |u| 统计：")
bypart = {}
for i, mo in enumerate(mag):
    bypart.setdefault(part_of.get(i + 1, "?"), []).append(mo)
for nm, v in sorted(bypart.items(), key=lambda kv: -max(kv[1])):
    v = np.array(v)
    print(f"    {nm:26s} n={len(v):>6d}  |u|max={v.max():.4f}  "
          f"中位={np.median(v):.4f}  p99={np.percentile(v,99):.4f}")
print(f"\n  支反力沿轴 = {(np.asarray(next((v for k,v in ds.items() if 'reaction' in k.lower()), 0))[-1].reshape(-1,3) @ ax).sum():.3f} N")
