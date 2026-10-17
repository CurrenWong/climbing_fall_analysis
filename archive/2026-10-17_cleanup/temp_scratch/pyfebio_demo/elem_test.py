"""单元级验证：单个楔形单向压缩，与解析解比。

底三角 (0,0,0)(10,0,0)(0,10,0)，顶面 z=5；固定底面，顶面压 0.1 mm。
解析：A = 50 mm², L = 5, ε = 0.1/5 = 0.02 ⇒ F = E·A·ε
      E = 10 MPa ⇒ F = 10 × 50 × 0.02 = **10 N**

对比：penta6 / penta15 / 3×tet4（WEDGE_T 拆法）/ hex8
若 penta6 明显低于 10 N ⇒ 该单元在 FEBio 里被欠积分（沙漏）
若 tet4 拆法明显高于 10 N ⇒ 过约束
"""
import subprocess
from pathlib import Path

TMP = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\elem_test")
TMP.mkdir(parents=True, exist_ok=True)
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
E = 10.0
NU = 0.0                      # ★ 用 ν=0 去掉体积耦合，纯看单元本身
H = 5.0
DZ = 0.1
A = 50.0
F_EXACT = E * A * (DZ / H)

BOT = [(0, 0, 0), (10, 0, 0), (0, 10, 0)]
TOP = [(0, 0, H), (10, 0, H), (0, 10, H)]


def mid(a, b):
    return tuple((a[i] + b[i]) / 2.0 for i in range(3))


VARIANTS = {
    "penta6": (6, [(0, 1, 2, 3, 4, 5)]),
    "tet4x3": (6, [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]),
    "tet4x6_hexstyle": (6, [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]),
    "penta15": (15, [(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)]),
}
ETYPE = {"penta6": "penta6", "tet4x3": "tet4",
         "tet4x6_hexstyle": "tet4", "penta15": "penta15"}

P15_MID = [mid(BOT[0], BOT[1]), mid(BOT[1], BOT[2]), mid(BOT[2], BOT[0]),
           mid(TOP[0], TOP[1]), mid(TOP[1], TOP[2]), mid(TOP[2], TOP[0]),
           mid(BOT[0], TOP[0]), mid(BOT[1], TOP[1]), mid(BOT[2], TOP[2])]


def build(nnode, elems, etype):
    pts = list(BOT) + list(TOP)
    if nnode == 15:
        pts += P15_MID
    nodes = "\n".join(f'      <node id="{i+1}">{p[0]},{p[1]},{p[2]}</node>'
                      for i, p in enumerate(pts))
    els = "\n".join(
        f'      <elem id="{i+1}">'
        f'{",".join(str(c+1) for c in e)}</elem>'
        for i, e in enumerate(elems))
    return f"""<?xml version="1.0" encoding="utf-8"?>
<febio_spec version="4.0">
  <Module type="solid"/>
  <Material>
    <material id="1" name="m" type="isotropic elastic">
      <density>1e-9</density><E>{E}</E><v>{NU}</v>
    </material>
  </Material>
  <Mesh>
    <Nodes name="all">
{nodes}
    </Nodes>
    <Elements type="{etype}" name="e">
{els}
    </Elements>
    <NodeSet name="bot">1,2,3</NodeSet>
    <NodeSet name="top">4,5,6</NodeSet>
  </Mesh>
  <MeshDomains>
    <SolidDomain name="e" mat="m"/>
  </MeshDomains>
  <Boundary>
    <bc name="fix" type="zero displacement" node_set="bot">
      <x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof>
    </bc>
    <bc name="push" type="prescribed displacement" node_set="top">
      <x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof>
      <value>0,0,-{DZ}</value>
    </bc>
  </Boundary>
  <Output>
    <plotfile type="febio"><var type="reaction forces"/></plotfile>
  </Output>
  <Control>
    <analysis>STATIC</analysis><time_steps>1</time_steps><step_size>1.0</step_size>
    <solver type="solid"><linear_solver type="pardiso"/></solver>
  </Control>
</febio_spec>
"""


print(f"解析解 F = E·A·ε = {E}×{A}×{DZ/H} = {F_EXACT:.3f} N   (ν={NU}, 无体积耦合)")
print(f"{'单元':>18} {'rc':>3} {'F(N)':>10} {'比值':>8}")
print("-" * 48)
for name, (nn, els) in VARIANTS.items():
    d = TMP / name
    d.mkdir(parents=True, exist_ok=True)
    (d / "t.feb").write_text(build(nn, els, ETYPE[name]), encoding="utf-8")
    r = subprocess.run([str(FEBIO), "-i", "t.feb"], cwd=d,
                       capture_output=True, text=True, timeout=180)
    out = (r.stdout or "") + (r.stderr or "")
    if "N O R M A L" not in out:
        err = [l.strip() for l in out.split("\n") if "ERROR" in l][:1]
        print(f"{name:>18} {'1':>3} {'-':>10} {'-':>8}  {err}")
        continue
    # 从 .log 或 xplt 读支反力
    import numpy as np
    try:
        import pyfebio.xplt as X
        xp = sorted(d.glob("*.xplt"))[0]
        X.to_hdf5(str(xp), str(d / "t.h5"))
        import h5py
        with h5py.File(d / "t.h5", "r") as h:
            st = sorted(int(k) for k in h["states"])
            R = np.asarray(
                h[f"states/{max(st)}/node_data/reaction forces/all"][()],
                float).reshape(-1, 3)
        Fz = abs(R[:, 2].sum())
        print(f"{name:>18} {0:>3} {Fz:10.3f} {Fz/F_EXACT:8.3f}")
    except Exception as e:
        print(f"{name:>18} {'?':>3} 读结果失败 {e}")
