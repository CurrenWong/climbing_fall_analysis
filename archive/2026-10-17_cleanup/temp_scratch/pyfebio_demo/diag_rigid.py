"""胫骨位移剖面：是弹性拉伸还是刚体模态？刚体拟合残差。"""
import inspect
import re
import textwrap
from pathlib import Path

import h5py
import numpy as np
import pyfebio.xplt as X

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio")
txt = (W / "_base.feb").read_text(encoding="latin-1")
parts = {}
for m in re.finditer(
        r'<Elements\b[^>]*type="(tet4|hex8|penta6)"[^>]*name="([^"]+)"[^>]*>(.*?)</Elements>',
        txt, re.S):
    parts[m.group(2)] = sorted({int(y) for s in
                                re.findall(r"<elem\b[^>]*>\s*([\d,]+)\s*</elem>", m.group(3))
                                for y in s.split(",")})
XYZ = {}
for m in re.finditer(
        r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)', txt):
    XYZ[int(m.group(1))] = np.array([float(m.group(2)), float(m.group(3)),
                                     float(m.group(4))])

for _o, _n in (('domain["name"]', 'domain.get("name", str(domain["id"][0]))'),
               ('surface["name"]', 'surface.get("name", str(surface["id"][0]))'),
               ("{domain['name']}", "{domain.get('name', str(domain['id'][0]))}"),
               ("{surface['name']}", "{surface.get('name', str(surface['id'][0]))}"),
               ('if nodeset["name"] == "":', 'if nodeset.get("name", "") == "":'),
               ("{nodeset['name']}", "{nodeset.get('name', str(nodeset['id'][0]))}"),
               ("{elementset['name']}",
                "{elementset.get('name', str(elementset['id'][0]))}")):
    pass
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
X.to_hdf5(str(W / "_base.xplt"), str(W / "_ba.h5"))
with h5py.File(W / "_ba.h5", "r") as h:
    ns = sorted((int(k) for k in h["states"]), key=int)
    U = np.asarray(h[f"states/{ns[-1]}/node_data/displacement/all"][()], float)
    # 找时间数据集
    print("states 数 =", len(ns))

AX = np.array([0.8153, 0.135, 0.5631])
AX = AX / np.linalg.norm(AX)

tb = parts["tibia_cort__hex8"]
P = np.array([XYZ[i] for i in tb])
D = U[[i - 1 for i in tb]]
proj = P @ AX

print("\n" + "=" * 74)
print("① 胫骨皮质：沿轴每 25mm 的位移（轴向分量 u·axis 与模 |u|）")
print("=" * 74)
print(f"  {'轴位置区间':>18} {'n':>6} {'u·axis 均值':>12} {'u·axis 范围':>24} {'|u| 均值':>10}")
lo, hi = proj.min(), proj.max()
for a in np.arange(np.floor(lo / 25) * 25, hi, 25):
    sel = (proj >= a) & (proj < a + 25)
    if sel.sum() < 5:
        continue
    ua = D[sel] @ AX
    print(f"  [{a:7.0f},{a+25:7.0f}) {sel.sum():>6d} {ua.mean():>12.4f} "
          f"[{ua.min():>10.4f},{ua.max():>10.4f}] "
          f"{np.linalg.norm(D[sel],axis=1).mean():>10.4f}")

print("\n" + "=" * 74)
print("② 刚体拟合：u ≈ a + ω×(x−x0)，看残差")
print("=" * 74)
X0 = P.mean(0)
# 构造 [I | -skew(x-x0)] 形：u = t + w × r  => u = t - skew(r) w
R = P - X0
M = np.zeros((3 * len(P), 6))
for k in range(len(P)):
    M[3 * k:3 * k + 3, 0:3] = np.eye(3)
    M[3 * k:3 * k + 3, 3:6] = -np.array([[0, -R[k, 2], R[k, 1]],
                                         [R[k, 2], 0, -R[k, 0]],
                                         [-R[k, 1], R[k, 0], 0]])
b = D.reshape(-1)
c, *_ = np.linalg.lstsq(M, b, rcond=None)
res = (M @ c).reshape(-1, 3) - D
print(f"  平移 t = {c[0:3]}")
print(f"  旋转 w = {c[3:6]}  (|w| = {np.linalg.norm(c[3:6]):.5f} rad "
      f"= {np.degrees(np.linalg.norm(c[3:6])):.4f}°)")
print(f"  刚体拟合残差: 中位 {np.linalg.norm(res,axis=1).mean():.5f}  "
      f"max {np.linalg.norm(res,axis=1).max():.5f} mm")
print(f"  原始 |u|:     中位 {np.linalg.norm(D,axis=1).mean():.5f}  "
      f"max {np.linalg.norm(D,axis=1).max():.5f} mm")
print(f"  ⇒ 残差/原始 = {np.linalg.norm(res,axis=1).max()/np.linalg.norm(D,axis=1).max()*100:.1f}%")

print("\n" + "=" * 74)
print("③ 相同刚体拟合对 —— 只取「整根胫骨」（含髓腔/松质）")
print("=" * 74)
alls = sorted({i for d in ("tibia_cort__hex8", "tibia_cort__tet4",
                           "tibia_end_spon__tet4", "tibia_marrow__tet4")
               for i in parts[d]})
Pa = np.array([XYZ[i] for i in alls])
Da = U[[i - 1 for i in alls]]
Ra = Pa - Pa.mean(0)
Ma = np.zeros((3 * len(Pa), 6))
for k in range(len(Pa)):
    Ma[3 * k:3 * k + 3, 0:3] = np.eye(3)
    Ma[3 * k:3 * k + 3, 3:6] = -np.array([[0, -Ra[k, 2], Ra[k, 1]],
                                          [Ra[k, 2], 0, -Ra[k, 0]],
                                          [-Ra[k, 1], Ra[k, 0], 0]])
ca, *_ = np.linalg.lstsq(Ma, Da.reshape(-1), rcond=None)
resa = (Ma @ ca).reshape(-1, 3) - Da
print(f"  n={len(Pa)}  残差 max {np.linalg.norm(resa,axis=1).max():.5f} mm "
      f"/ 原始 max {np.linalg.norm(Da,axis=1).max():.5f} mm "
      f"= {np.linalg.norm(resa,axis=1).max()/np.linalg.norm(Da,axis=1).max()*100:.1f}%")
print(f"  旋转 |w| = {np.linalg.norm(ca[3:6]):.5f} rad = "
      f"{np.degrees(np.linalg.norm(ca[3:6])):.4f}°")
