"""探 FEBio 4.13 是否支持 penta15，以及其节点序。

VTK 二次楔形序（待验证）：
  0,1,2 底三角 | 3,4,5 顶三角
  6=mid(0,1) 7=mid(1,2) 8=mid(2,0)
  9=mid(3,4) 10=mid(4,5) 11=mid(5,3)
  12=mid(0,3) 13=mid(1,4) 14=mid(2,5)
"""
import subprocess
from pathlib import Path

TMP = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\p15_probe")
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")

C = [(0, 0, 0), (2, 0, 0), (0, 2, 0),          # 底三角
     (0, 0, 1), (2, 0, 1), (0, 2, 1)]          # 顶三角
V = C + []                                      # 15 个节点：6 角 + 9 边中点


def mid(a, b):
    return tuple((C[a][i] + C[b][i]) / 2.0 for i in range(3))


ORDERS = {
    "VTK": C + [mid(0, 1), mid(1, 2), mid(2, 0),
                mid(3, 4), mid(4, 5), mid(5, 3),
                mid(0, 3), mid(1, 4), mid(2, 5)],
    "ALT_bottom_first_edges": C + [mid(0, 1), mid(1, 2), mid(2, 0),
                                   mid(0, 3), mid(1, 4), mid(2, 5),
                                   mid(3, 4), mid(4, 5), mid(5, 3)],
}


def build(nodes):
    n = "\n".join(f'      <node id="{i+1}">{p[0]},{p[1]},{p[2]}</node>'
                  for i, p in enumerate(nodes))
    return f"""<?xml version="1.0" encoding="utf-8"?>
<febio_spec version="4.0">
  <Module type="solid"/>
  <Material>
    <material id="1" name="m" type="isotropic elastic">
      <density>1e-9</density><E>10.0</E><v>0.3</v>
    </material>
  </Material>
  <Mesh>
    <Nodes name="all">
{n}
    </Nodes>
    <Elements type="penta15" name="e">
      <elem id="1">1,2,3,4,5,6,7,8,9,10,11,12,13,14,15</elem>
    </Elements>
    <NodeSet name="n1">1,2,3,4,5,6,7,8,9,10,11,12,13,14,15</NodeSet>
  </Mesh>
  <MeshDomains>
    <SolidDomain name="e" mat="m"/>
  </MeshDomains>
  <Boundary>
    <bc name="fix" type="zero displacement" node_set="n1">
      <x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof>
    </bc>
  </Boundary>
  <Control>
    <analysis>STATIC</analysis><time_steps>1</time_steps><step_size>1.0</step_size>
    <solver type="solid"><linear_solver type="pardiso"/></solver>
  </Control>
</febio_spec>
"""


for name, nodes in ORDERS.items():
    d = TMP / name
    d.mkdir(parents=True, exist_ok=True)
    (d / "t.feb").write_text(build(nodes), encoding="utf-8")
    r = subprocess.run([str(FEBIO), "-i", "t.feb"], cwd=d,
                       capture_output=True, text=True, timeout=180)
    out = (r.stdout or "") + (r.stderr or "")
    errs = [l.strip() for l in out.split("\n")
            if ("ERROR" in l or "unrecognized" in l or "Unknown" in l
                or "invalid" in l.lower()
                or "Negative jacobian" in l or "FAILED" in l)]
    ok = "N O R M A L" in out
    print(f"  {name:26s} {'OK 正常终止' if ok else '? ' + str(errs[:3])}")
