"""最小实验（修正版）：确定 FEBio 4.13 的 penta6 节点序约定。"""
import subprocess
from pathlib import Path

TMP = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\penta_probe")
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")

COORD = [(0, 0, 0), (1, 0, 0), (0, 1, 0),
         (0, 0, 1), (1, 0, 1), (0, 1, 1)]

VARIANTS = {
    "A": [1, 2, 3, 4, 5, 6],
    "B": [1, 2, 3, 6, 5, 4],
    "C": [1, 2, 3, 5, 6, 4],
    "D": [1, 2, 3, 4, 6, 5],
    "E": [1, 3, 2, 4, 5, 6],
    "F": [2, 1, 3, 4, 5, 6],
}


def build(conn):
    nodes = "\n".join(f'      <node id="{i+1}">{x},{y},{z}</node>'
                      for i, (x, y, z) in enumerate(COORD))
    return f"""<?xml version="1.0" encoding="utf-8"?>
<febio_spec version="4.0">
  <Module type="solid"/>
  <Material>
    <material id="1" name="m" type="isotropic elastic">
      <density>1e-9</density>
      <E>10.0</E>
      <v>0.4</v>
    </material>
  </Material>
  <Mesh>
    <Nodes name="all">
{nodes}
    </Nodes>
    <Elements type="penta6" name="e">
      <elem id="1">{",".join(str(c) for c in conn)}</elem>
    </Elements>
    <NodeSet name="n1">1,2,3,4,5,6</NodeSet>
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
    <analysis>STATIC</analysis>
    <time_steps>1</time_steps>
    <step_size>1.0</step_size>
    <solver type="solid">
      <linear_solver type="pardiso"/>
    </solver>
  </Control>
</febio_spec>
"""


for name, conn in VARIANTS.items():
    d = TMP / name
    d.mkdir(parents=True, exist_ok=True)
    (d / "t.feb").write_text(build(conn), encoding="utf-8")
    r = subprocess.run([str(FEBIO), "-i", "t.feb"], cwd=d,
                       capture_output=True, text=True, timeout=180)
    out = (r.stdout or "") + (r.stderr or "")
    if "Negative jacobian" in out:
        verdict = "X 负 Jacobian(初始化)"
    elif "FAILED!" in out and "Reading file" in out:
        err = [l for l in out.split("\n") if "unrecognized" in l
               or "tag" in l.lower()][:1]
        verdict = f"? parse: {err}"
    elif "N O R M A L" in out:
        verdict = "OK 正常终止"
    elif "ERROR" in out:
        err = [l.strip() for l in out.split("\n")
               if "ERROR" in l or "error" in l][:2]
        verdict = f"? {err}"
    else:
        verdict = "? 未判"
    print(f"  {name} {str(conn):26s} {verdict}")
