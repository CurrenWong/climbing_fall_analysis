"""truss 写法矩阵探测：BeamDomain type × 材料 × 面积写法。"""
import subprocess
from pathlib import Path

WORK = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\truss_probe")
WORK.mkdir(parents=True, exist_ok=True)
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")

HEAD = '''<?xml version="1.0" encoding="ISO-8859-1"?>
<febio_spec version="4.0">
  <Module type="solid"/>
  <Globals><Constants><T>0</T><R>0</R><Fc>0</Fc></Constants></Globals>
  <Material>
{mat}
  </Material>
  <Mesh>
    <Nodes name="all">
      <node id="1">0,0,0</node>
      <node id="2">0,0,-10</node>
    </Nodes>
    <Elements type="line2" name="tr">
      <elem id="1">1,2</elem>
    </Elements>
    <NodeSet name="fix">1</NodeSet>
    <NodeSet name="pull">2</NodeSet>
  </Mesh>
  <MeshDomains>
    <BeamDomain name="tr" mat="lig" type="{dt}"{extra}/>
  </MeshDomains>
  <LoadData>
    <load_controller id="1" type="loadcurve">
      <interpolate>LINEAR</interpolate><extend>CONSTANT</extend>
      <points><pt>0.0,0.0</pt><pt>1.0,1.0</pt></points>
    </load_controller>
  </LoadData>
  <Boundary>
    <bc name="f" type="zero displacement" node_set="fix">
      <x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof>
    </bc>
    <bc name="p" type="prescribed displacement" node_set="pull">
      <dof>z</dof><value lc="1">-1.0</value><relative>0</relative>
    </bc>
  </Boundary>
  <Output>
    <plotfile type="febio"><var type="displacement"/><var type="reaction forces"/></plotfile>
  </Output>
  <Control>
    <analysis>STATIC</analysis><time_steps>2</time_steps><step_size>0.5</step_size>
    <solver type="solid"><linear_solver type="pardiso"/></solver>
  </Control>
</febio_spec>
'''

M_ISO = '''    <material id="1" name="lig" type="isotropic elastic">
      <density>1</density><E>260</E><v>0.3</v>
    </material>'''
M_ISO_AREA = '''    <material id="1" name="lig" type="isotropic elastic">
      <density>1</density><E>260</E><v>0.3</v>
      <area>10</area>
    </material>'''
M_LE = '''    <material id="1" name="lig" type="linear elastic">
      <density>1</density><E>260</E><v>0.3</v>
    </material>'''
M_RIGID = '''    <material id="1" name="lig" type="rigid">
      <density>1</density>
    </material>'''

CASES = [
    ("et_iso", "elastic-truss", M_ISO, ""),
    ("et_iso_area", "elastic-truss", M_ISO_AREA, ""),
    ("lt_le", "linear-truss", M_LE, ""),
    ("lt_rigid", "linear-truss", M_RIGID, ""),
    ("lb_iso", "linear-beam", M_ISO, ""),
]

for name, dt, mat, extra in CASES:
    f = WORK / f"{name}.feb"
    f.write_text(HEAD.format(dt=dt, mat=mat, extra=extra), encoding="latin-1")
    r = subprocess.run([str(FEBIO), "-i", str(f)], cwd=str(WORK),
                       capture_output=True, text=True, timeout=300)
    o = (r.stdout or "") + (r.stderr or "")
    ok = "NORMAL TERMINATION" in o
    print(f"\n--- {name} (BeamDomain type={dt}) rc={r.returncode} "
          f"{'** NORMAL TERMINATION **' if ok else 'FAILED'} ---")
    tail = [l.strip() for l in o.split("\n") if l.strip()][-6:]
    for l in tail:
        print("   ", l[:130])
