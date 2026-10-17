"""A/B 测：FEBio 4.13 里 surface 集合该用 <Elements type="quad4"> 还是 <Surface>？"""
import subprocess
from pathlib import Path

FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio\probe")
W.mkdir(parents=True, exist_ok=True)

HEAD = """<?xml version="1.0" encoding="ISO-8859-1"?>
<febio_spec version="4.0">
  <Module type="solid"/>
  <Globals><Constants><T>298.0</T><R>8.314e-06</R><Fc>9.6485e-05</Fc></Constants></Globals>
  <Material>
    <material name="m" id="1" type="isotropic elastic">
      <density>1e-06</density><E>1000</E><v>0.3</v>
    </material>
  </Material>
  <Mesh>
    <Nodes name="nodes">
      <node id="1">0,0,0</node><node id="2">1,0,0</node>
      <node id="3">1,1,0</node><node id="4">0,1,0</node>
      <node id="5">0,0,1</node><node id="6">1,0,1</node>
      <node id="7">1,1,1</node><node id="8">0,1,1</node>
    </Nodes>
    <Elements name="box" type="hex8"><elem id="1">1,2,3,4,5,6,7,8</elem></Elements>
"""
TAILSURF = {
    "A_Elements_quad4": '    <Elements name="topsurf" type="quad4"><elem id="1">5,6,7,8</elem></Elements>\n',
    "B_Surface_tag": '    <Surface name="topsurf"><elem id="1">5,6,7,8</elem></Surface>\n',
    "C_Surface_with_type": '    <Surface name="topsurf" type="quad4"><elem id="1">5,6,7,8</elem></Surface>\n',
}
REST = """    <NodeSet name="bottom">1,2,3,4</NodeSet>
  </Mesh>
  <MeshDomains><SolidDomain name="box" mat="m"/></MeshDomains>
  <Loads>
    <surface_load type="traction" surface="topsurf">
      <traction>0,0,-1</traction>
    </surface_load>
  </Loads>
  <Boundary>
    <bc type="zero displacement" node_set="bottom">
      <x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof>
    </bc>
  </Boundary>
  <Control>
    <analysis>STATIC</analysis><time_steps>2</time_steps><step_size>0.5</step_size>
    <solver type="solid"><linear_solver type="pardiso"/></solver>
  </Control>
</febio_spec>
"""

for name, blk in TAILSURF.items():
    f = W / f"{name}.feb"
    f.write_text(HEAD + blk + REST, encoding="latin-1")
    r = subprocess.run([str(FEBIO), "-i", str(f)], cwd=str(W),
                       capture_output=True, text=True, timeout=300)
    out = (r.stdout or "") + (r.stderr or "")
    err = ""
    if "ERROR" in out:
        i = out.find("ERROR")
        err = " | ".join(x.strip(" *") for x in out[i:i+260].split("\n")[1:4] if x.strip())
    print(f"{'OK  ' if r.returncode == 0 else 'FAIL'} {name:22s} rc={r.returncode}  {err[:120]}")
    if r.returncode == 0:
        print("       ->", [l.strip() for l in out.split("\n") if "NORMAL" in l or "error" in l.lower()][:3])
