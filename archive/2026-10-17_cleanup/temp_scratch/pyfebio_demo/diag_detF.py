"""诊断：逐单元算 det(F)，找濒临/已经翻转的单元（定位到域）。

原理：FEBio 的 "Negative jacobian" = 某单元变形梯度 det(F) <= 0。
末个收敛态的位移已知 ⇒ 可算出每个单元的 det(F)，最小值即下一个翻转者。
"""
import re
import sys
import textwrap
import inspect
import h5py
import numpy as np
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
W = ROOT / "temp/thums/febio_lig"
FEB = W / "ankle_lig_WK.feb"

t = FEB.read_text(encoding="utf-8")

# ---- 节点 ----
nsec = re.search(r"<Nodes[^>]*>(.*?)</Nodes>", t, re.S).group(1)
nid, coords = [], []
for m in re.finditer(r'<node id="(\d+)">([^<]+)</node>', nsec):
    v = [float(x) for x in re.split(r"[,\s]+", m.group(2).strip()) if x]
    nid.append(int(m.group(1)))
    coords.append(v[:3])
nid = np.array(nid)
X0 = np.array(coords)
lut = {n: k for k, n in enumerate(nid)}
print(f"节点 {len(nid)}")

# ---- 单元（按 SolidDomain 分组）----
esec = re.search(r"<Elements>(.*?)</Elements>", t, re.S)
if not esec:
    print("!! 找不到 <Elements>")
    print([m.group(0) for m in re.finditer(r"<[A-Za-z_]+>", t)][:40])
    sys.exit(1)
EB = esec.group(1)
print("元素标签：", sorted({m.group(1) for m in re.finditer(
    r"<(tet4|hex8|tet10|hex20|penta6|tri3|quad4|beam2|truss2)[ >]", EB)}))

DOM = {}
# 按行扫描：<SolidDomain name=.../> 之后到下一个 domain 前的元素归它
cur = None
for line in EB.split("\n"):
    m = re.search(r'<SolidDomain name="([^"]+)"', line)
    if m:
        cur = m.group(1)
        DOM.setdefault(cur, [])
        continue
    m = re.search(r"<(tet4|hex8|tet10|penta6)[^>]*>([^<]+)</", line)
    if m and cur:
        ids = [int(x) for x in re.split(r"[,\s]+", m.group(2).strip()) if x]
        DOM[cur].append((m.group(1), ids))
print(f"\n域 {len(DOM)} 个")
for k, v in DOM.items():
    print(f"  {k:34s} {len(v):6d} 单元 ({v[0][0] if v else '-'})")

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

src = sorted(W.glob("ankle_lig_WK*.xplt"),
             key=lambda p: p.stat().st_mtime, reverse=True)[0]
dst = W / "_dg.h5"
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
D = np.asarray(disp, float)
if D.ndim == 3:
    D = D[-1]
D = D.reshape(-1, 3)
_c = np.asarray(coords_plt)
Xp = (np.column_stack([_c["x"], _c["y"], _c["z"]]).astype(float)
      if _c.dtype.names else _c.reshape(-1, 3).astype(float))
Xd = Xp + D
print(f"\n末帧 |u|max = {np.linalg.norm(D, axis=1).max():.4f} mm")

# ---- det(F) ----
def detF_tet(idx):
    """tet4: 用 3 条边做变形梯度。"""
    x = X0[idx]                    # (N,4,3)
    y = Xd[idx]
    F = np.stack([y[:, 1] - y[:, 0], y[:, 2] - y[:, 0], y[:, 3] - y[:, 0]], 1)
    G = np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]], 1)
    return np.linalg.det(F) / np.linalg.det(G)


def detF_hex(idx):
    """hex8: 8 个角点处的 det(J)（Gauss 点近似用角点）。"""
    # 在中心点算（用对面的差）
    x = X0[idx]
    y = Xd[idx]
    e = [(1, 0), (2, 3), (4, 7), (5, 6)]     # 近似三对边
    F = np.stack([y[:, e[0][1]] - y[:, e[0][0]],
                  y[:, e[1][1]] - y[:, e[1][0]],
                  y[:, e[2][1]] - y[:, e[2][0]]], 1)
    G = np.stack([x[:, e[0][1]] - x[:, e[0][0]],
                  x[:, e[1][1]] - x[:, e[1][0]],
                  x[:, e[2][1]] - x[:, e[2][0]]], 1)
    return np.linalg.det(F) / np.linalg.det(G)


print("\n" + "=" * 82)
print(f"{'域':34s} {'单元':>6} {'det(F) min':>12} {'p1':>10} {'p5':>10} {'<0 个数':>8}")
print("-" * 82)
worst = []
for k, els in DOM.items():
    et = els[0][0]
    idx = np.array([[lut[n] for n in ids] for _, ids in els])
    d = detF_tet(idx) if et == "tet4" else detF_hex(idx)
    worst.append((d.min(), k, et, int((d <= 0).sum()), len(d)))
    print(f"{k:34s} {len(d):6d} {d.min():12.5f} "
          f"{np.percentile(d, 1):10.5f} {np.percentile(d, 5):10.5f} "
          f"{int((d <= 0).sum()):8d}")
worst.sort()
print("\n最危险的 6 个域：")
for mn, k, et, nneg, n in worst[:6]:
    print(f"  {k:34s} ({et}) det min {mn:10.5f}  负 {nneg}/{n}")
