"""单元级验证（修正 BC 语法）：单个楔形单向压缩 vs 解析解。

★ FEBio 4.13 的 prescribed displacement 语法（从主模型抄来的已验证写法）：
    <bc type="prescribed displacement" node_set="...">
      <dof>x</dof><value>-0.1</value><relative>0</relative>
    </bc>
  每个自由度一条 BC；**不是** <x_dof>1</x_dof>（那是 zero displacement 的写法）。
"""
import subprocess
from pathlib import Path
import numpy as np
import h5py

TMP = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\elem_test2")
TMP.mkdir(parents=True, exist_ok=True)
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
E, NU, H, DZ, A = 10.0, 0.0, 5.0, 0.1, 50.0
F_EXACT = E * A * (DZ / H)

BOT = [(0, 0, 0), (10, 0, 0), (0, 10, 0)]
TOP = [(0, 0, H), (10, 0, H), (0, 10, H)]


def mid(a, b):
    return tuple((a[i] + b[i]) / 2.0 for i in range(3))


P15 = [mid(BOT[0], BOT[1]), mid(BOT[1], BOT[2]), mid(BOT[2], BOT[0]),
       mid(TOP[0], TOP[1]), mid(TOP[1], TOP[2]), mid(TOP[2], TOP[0]),
       mid(BOT[0], TOP[0]), mid(BOT[1], TOP[1]), mid(BOT[2], TOP[2])]

# (名称, 节点数, FEBio 单元类型, 连接表（0-based）)
VARIANTS = [
    ("penta6", 6, "penta6", [(0, 1, 2, 3, 4, 5)]),
    ("tet4x3", 6, "tet4", [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]),
    ("penta15", 15, "penta15",
     [(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)]),
]


def build(nnode, etype, elems):
    pts = list(BOT) + list(TOP) + (P15 if nnode == 15 else [])
    nodes = "\n".join(f'      <node id="{i+1}">{p[0]},{p[1]},{p[2]}</node>'
                      for i, p in enumerate(pts))
    els = "\n".join(f'      <elem id="{i+1}">'
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
    <bc name="pushx" type="prescribed displacement" node_set="top">
      <dof>x</dof><value>0</value><relative>0</relative>
    </bc>
    <bc name="pushy" type="prescribed displacement" node_set="top">
      <dof>y</dof><value>0</value><relative>0</relative>
    </bc>
    <bc name="pushz" type="prescribed displacement" node_set="top">
      <dof>z</dof><value>-{DZ}</value><relative>0</relative>
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


print(f"解析解 F = E·A·ε = {E}×{A}×{DZ/H} = {F_EXACT:.3f} N  (ν={NU}，无体积耦合)")
print(f"{'单元':>10} {'rc':>3} {'F_z(N)':>10} {'F/F_exact':>10}  备注")
print("-" * 56)
import pyfebio.xplt as X
for name, nn, et, els in VARIANTS:
    d = TMP / name
    d.mkdir(parents=True, exist_ok=True)
    (d / "t.feb").write_text(build(nn, et, els), encoding="utf-8")
    r = subprocess.run([str(FEBIO), "-i", "t.feb"], cwd=d,
                       capture_output=True, text=True, timeout=180)
    o = (r.stdout or "") + (r.stderr or "")
    if "N O R M A L" not in o:
        e = [l.strip() for l in o.split("\n")
             if "ERROR" in l or "unrecognized" in l][:1]
        print(f"{name:>10} {'1':>3} {'-':>10} {'-':>10}  {e}")
        continue
    xp = sorted(d.glob("*.xplt"))[0]
    X.to_hdf5(str(xp), str(d / "t.h5"))
    with h5py.File(d / "t.h5", "r") as h:
        st = sorted(int(k) for k in h["states"])
        R = np.asarray(
            h[f"states/{max(st)}/node_data/reaction forces/all"][()],
            float).reshape(-1, 3)
    Fz = abs(R[:, 2].sum())
    print(f"{name:>10} {0:>3} {Fz:10.3f} {Fz/F_EXACT:10.3f}")
