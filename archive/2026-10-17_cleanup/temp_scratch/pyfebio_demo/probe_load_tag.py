"""穷举 FEBio 4.13 的 load 段写法，用最小立方体模型让 FEBio 自己判断。
基线用已知可跑通的 temp/pyfebio_demo/my_model.feb。
"""
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
    <Elements name="topsurf" type="quad4"><elem id="1">5,6,7,8</elem></Elements>
    <NodeSet name="bottom">1,2,3,4</NodeSet>
    <NodeSet name="topsurf_ns">5,6,7,8</NodeSet>
  </Mesh>
  <MeshDomains><SolidDomain name="box" mat="m"/></MeshDomains>
  <Boundary>
    <bc type="zero displacement" node_set="bottom">
      <x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof>
    </bc>
  </Boundary>
"""

TAIL = """  <Control>
    <analysis>STATIC</analysis><time_steps>3</time_steps><step_size>0.33</step_size>
    <solver type="solid"><linear_solver type="pardiso"/></solver>
  </Control>
</febio_spec>
"""

VARIANTS = {
    "v1_Loads_load_surface_traction": """  <Loads>
    <load name="L" type="surface traction">
      <surface>topsurf</surface><traction>0,0,-1</traction>
    </load>
  </Loads>
""",
    "v2_Loads_load_surface_traction_underscore": """  <Loads>
    <load name="L" type="surface_traction">
      <surface>topsurf</surface><traction>0,0,-1</traction>
    </load>
  </Loads>
""",
    "v3_lowercase_loads": """  <loads>
    <load name="L" type="surface traction">
      <surface>topsurf</surface><traction>0,0,-1</traction>
    </load>
  </loads>
""",
    "v4_Loads_surface_load": """  <Loads>
    <surface_load name="L" type="surface traction">
      <surface>topsurf</surface><traction>0,0,-1</traction>
    </surface_load>
  </Loads>
""",
    "v5_Loads_nodal_force": """  <Loads>
    <load name="L" type="nodal force">
      <node_set>topsurf_ns</node_set><force>0,0,-1</force>
    </load>
  </Loads>
""",
    "v6_no_bc_only_pressure": """  <Loads>
    <load name="L" type="pressure">
      <surface>topsurf</surface><pressure>1.0</pressure>
    </load>
  </Loads>
""",
}

for name, block in VARIANTS.items():
    f = W / f"{name}.feb"
    f.write_text(HEAD + block + TAIL, encoding="latin-1")
    r = subprocess.run([str(FEBIO), "-i", str(f)], cwd=str(W),
                       capture_output=True, text=True, timeout=300)
    out = (r.stdout or "") + (r.stderr or "")
    err = ""
    if "ERROR" in out:
        i = out.find("ERROR")
        err = " | ".join(x.strip(" *") for x in out[i:i+300].split("\n")[1:4] if x.strip())
    ok = r.returncode == 0
    print(f"{'OK ' if ok else 'FAIL'} {name:44s} rc={r.returncode}  {err[:110]}")
