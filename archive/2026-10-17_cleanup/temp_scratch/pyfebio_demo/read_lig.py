"""读踝韧带模型的求解结果：每条韧带的位移（拉伸长度）+ 支反力。

对标吴恺的「韧带位移（拉伸长度）」：
  对每条韧带的每根"纤维"(a,b)，ΔL = |x_b+u_b - (x_a+u_a)| - |x_b - x_a|
  韧带的拉伸 = 各纤维 ΔL 的平均（并报最大）。
"""
import sys
import h5py
import numpy as np
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
W = ROOT / "temp/thums/febio_lig"

tag = sys.argv[1] if len(sys.argv) > 1 else "WK"
r = np.load(W / f"lig_read_{tag}.npz", allow_pickle=True)
feb = W / f"ankle_lig_{tag}.feb"

# 取 .xplt（最新）
xplts = sorted(W.glob(f"ankle_lig_{tag}*.xplt"),
               key=lambda p: p.stat().st_mtime, reverse=True)
if not xplts:
    print("找不到 .xplt")
    sys.exit(1)
src = xplts[0]
dst = W / "lig_result.h5"
print(f"读 {src.name} ({src.stat().st_size/1024/1024:.1f} MB)")

import textwrap
import inspect
import pyfebio.xplt as X

# ── 容错补丁：pyfebio 的 parse_mesh 对 **shell/beam 域**没有 "name" 字段，
#   会 KeyError 'name'（MEMORY 里记为「需容错补丁」）。用 id 兜底。──
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
exec(compile(_src, "<patched_parse_mesh>", "exec"), _ns)
X.parse_mesh = _ns["parse_mesh"]

# ── 第二个容错补丁：parse_state 对**自定义 Surface 集合**会 KeyError
#   （MEMORY 已记录该坑）—— 找不到 set_id 时用占位名兜底。──
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
if 'mesh_dict[set_lut[key]].get(set_id)' not in _ps:
    print("!! parse_state 补丁未命中，需人工检查")
_ns2 = X.__dict__
exec(compile(_ps, "<patched_parse_state>", "exec"), _ns2)
X.parse_state = _ns2["parse_state"]
print("pyfebio parse_state 已打容错补丁")
print("pyfebio parse_mesh 已打容错补丁")

X.to_hdf5(str(src), str(dst))

coords = None
disp = None
react = None
with h5py.File(dst, "r") as h:
    keys = []

    def walk(g, pre=""):
        for k in g:
            if isinstance(g[k], h5py.Group):
                walk(g[k], pre + "/" + k)
            else:
                keys.append(pre + "/" + k)
    walk(h)
    nk = [k for k in keys if k.endswith("/nodes") or "nodes" in k.lower()]
    for k in keys:
        kl = k.lower()
        if kl.endswith("/nodes"):
            coords = h[k][()]
        elif "displacement" in kl:
            disp = h[k][()]
        elif "reaction" in kl:
            react = h[k][()]
print(f"节点 {None if coords is None else coords.shape}  "
      f"位移 {None if disp is None else disp.shape}  "
      f"支反力 {None if react is None else react.shape}")
if disp is None:
    print("无位移数据；keys 前 30:", keys[:30])
    sys.exit(1)

# disp: (n_steps, n_nodes, 3?)；取末帧
D = np.asarray(disp, dtype=float)
if D.ndim == 3:
    D = D[-1]                       # 末帧
D = D.reshape(-1, 3)
print(f"位移矩阵形状 {D.shape}")

_c = np.asarray(coords)
if _c.dtype.names:                       # 结构化数组 (id, x, y, z)
    X0 = np.column_stack([_c["x"], _c["y"], _c["z"]]).astype(float)
else:
    X0 = _c.reshape(-1, 3).astype(float)  # (n, 3)
assert D.shape[0] == X0.shape[0], f"节点数不符 {D.shape} vs {X0.shape}"

pairs = r["lig_pairs"]
names = [str(x) for x in r["lig_names"]]
eas = r["lig_ea"]
npair = r["lig_npair"]

print("\n" + "=" * 78)
print("每条韧带的位移（拉伸长度，mm）")
print("=" * 78)
print(f"{'韧带':34s} {'纤维':>4} {'平均ΔL':>9} {'最大ΔL':>9} "
      f"{'EA':>7} {'|u|均':>9} {'|u|max':>9}")
print("-" * 96)
res = []
for i, nm in enumerate(names):
    ds = []
    for a, b in pairs[i][:npair[i]]:
        a, b = int(a) - 1, int(b) - 1
        L0 = np.linalg.norm(X0[b] - X0[a])
        L1 = np.linalg.norm((X0[b] + D[b]) - (X0[a] + D[a]))
        ds.append(L1 - L0)
    if not ds:
        continue
    ds = np.array(ds)
    # 口径二：该韧带**全部节点的位移模**（吴恺"位移"可能指这个）
    nd = sorted({int(x) - 1 for pr in pairs[i][:npair[i]] for x in pr})
    umag = float(np.linalg.norm(D[nd], axis=1).mean()) if nd else 0.0
    umax = float(np.linalg.norm(D[nd], axis=1).max()) if nd else 0.0
    res.append((float(ds.mean()), nm, len(ds), float(ds.max()),
                float(eas[i]), umag, umax))
    print(f"{nm:34s} {len(ds):>4} {ds.mean():>9.4f} {ds.max():>9.4f} "
          f"{eas[i]:>8.1f} {umag:>9.4f} {umax:>9.4f}")

print("\n" + "=" * 78)
print("排序（按平均拉伸位移，降序）")
print("=" * 78)
res.sort(reverse=True)
print("\n【口径一 ΔL = 拉伸长度变化】")
for k, (m, nm, n, mx, ea, ug, ux) in enumerate(res, 1):
    print(f"  {k:2d}. {nm:34s} {m:8.4f} mm   (EA={ea:.1f})")
print("\n【口径二 |u| = 韧带节点绝对位移均值】")
for k, (_, nm, _, _, ea, ug, ux) in enumerate(
        sorted(res, key=lambda t: -t[5]), 1):
    print(f"  {k:2d}. {nm:34s} {ug:8.4f} mm   (EA={ea:.1f})")

WU = {"MD2_tibiocalcaneal": "胫跟", "MD3_tibionavicular": "胫舟",
      "CFL_calcaneofibular": "跟腓", "SYND_tibiofibular": "下胫腓前"}
sub = [(m, WU[nm]) for m, nm, *_ in res if nm in WU]
usub = [(ug, WU[nm]) for _, nm, _, _, _, ug, _ in sorted(
    res, key=lambda t: -t[5]) if nm in WU]
if sub:
    print("\n吴恺要的四条（吴恺: 胫跟 5.24 > 胫舟 5.04 > 跟腓 2.05 > 下胫腓前 1.92）")
    print("  口径一 ΔL    : " + " > ".join(nm for _, nm in sub))
    print("  口径二 |u|   : " + " > ".join(nm for _, nm in usub))
    for m, nm in sub:
        print(f"    {nm:8s} ΔL={m:8.4f} mm")

if react is not None:
    R = np.asarray(react, dtype=float).reshape(-1, 3)
    kax = int(np.argmax(np.abs(r["axis"])))
    print(f"\n支反力沿加载轴分量之和 = {R[:, kax].sum():.3f}")
