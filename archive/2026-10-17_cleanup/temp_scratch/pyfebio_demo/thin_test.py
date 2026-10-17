"""薄楔形单元测试：匹配软骨的真实长宽比（厚 ~0.7 mm / 面内 ~4 mm）。

若薄单元下 penta6 与 tet4x3 开始分歧 ⇒ 差异源于单元格式在**扁平/歪斜**时的行为。
同时扫「歪斜度」：把顶面平移，模拟偏移方向与面不垂直。
"""
import subprocess
from pathlib import Path
import numpy as np
import h5py
import pyfebio.xplt as X

TMP = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\thin_test")
TMP.mkdir(parents=True, exist_ok=True)
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
E, NU = 10.0, 0.0


def build(h, side, skew, etype, elems, nn=6, extra=None):
    BOT = [(0, 0, 0), (side, 0, 0), (0, side, 0)]
    TOP = [(0 + skew[0], skew[1], h), (side + skew[0], skew[1], h),
           (0 + skew[0], side + skew[1], h)]
    pts = BOT + TOP + (extra or [])
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
    <bc name="px" type="prescribed displacement" node_set="top">
      <dof>x</dof><value>{skew[0]}</value><relative>0</relative>
    </bc>
    <bc name="py" type="prescribed displacement" node_set="top">
      <dof>y</dof><value>{skew[1]}</value><relative>0</relative>
    </bc>
    <bc name="pz" type="prescribed displacement" node_set="top">
      <dof>z</dof><value>{-dz}</value><relative>0</relative>
    </bc>
  </Boundary>
  <Output><plotfile type="febio"><var type="reaction forces"/></plotfile></Output>
  <Control>
    <analysis>STATIC</analysis><time_steps>1</time_steps><step_size>1.0</step_size>
    <solver type="solid"><linear_solver type="pardiso"/></solver>
  </Control>
</febio_spec>
"""


TET3 = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]

print("薄楔形：厚 h，底三角边长 side，顶面平移 skew")
print(f"{'h':>6} {'side':>6} {'skew':>12} {'应变':>7} "
      f"{'penta6 F':>10} {'tet4x3 F':>10} {'p/t':>7}")
print("-" * 74)
for h, side, skew in ((5.0, 10.0, (0, 0)),
                      (0.7, 4.0, (0, 0)),
                      (0.7, 4.0, (0.6, 0.0)),
                      (0.7, 4.0, (0.9, 0.4)),
                      (0.2, 4.0, (1.5, 0.5))):
    dz = 0.1 * h                       # 10% 轴向应变
    A = 0.5 * side * side
    Fexact = E * A * (dz / h)
    res = {}
    for nm, et, els in (("penta6", "penta6", [(0, 1, 2, 3, 4, 5)]),
                        ("tet4x3", "tet4", TET3)):
        d = TMP / f"{nm}_{h}_{side}_{skew[0]}_{skew[1]}"
        d.mkdir(parents=True, exist_ok=True)
        (d / "t.feb").write_text(build(h, side, skew, et, els),
                                 encoding="utf-8")
        r = subprocess.run([str(FEBIO), "-i", "t.feb"], cwd=d,
                           capture_output=True, text=True, timeout=180)
        o = (r.stdout or "") + (r.stderr or "")
        if "N O R M A L" not in o:
            res[nm] = float("nan")
            continue
        xp = sorted(d.glob("*.xplt"))[0]
        X.to_hdf5(str(xp), str(d / "t.h5"))
        with h5py.File(d / "t.h5", "r") as hh:
            st = sorted(int(k) for k in hh["states"])
            R = np.asarray(
                hh[f"states/{max(st)}/node_data/reaction forces/all"][()],
                float).reshape(-1, 3)
        res[nm] = abs(R[:, 2].sum())
    p, t = res.get("penta6", float("nan")), res.get("tet4x3", float("nan"))
    rt = p / t if (t and t == t) else float("nan")
    print(f"{h:6.2f} {side:6.1f} {str(skew):>12} {0.10:7.2f} "
          f"{p:10.3f} {t:10.3f} {rt:7.3f}   (解析 {Fexact:.3f})")
