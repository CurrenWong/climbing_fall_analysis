"""最小实验：摸清 FEBio truss 的 .feb 写法（BeamDomain + 材料/面积）。"""
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
    <material id="1" name="lig" type="{mt}">
{params}
    </material>
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
    <BeamDomain name="tr" mat="lig" type="linear-truss"{extra}/>
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

VARIANTS = {
    "A_iso": ("isotropic elastic",
                         "      <density>1</density><E>260</E><v>0.3</v>", ""),
    "B_iso_area": ("isotropic elastic",
                              "      <density>1</density><E>260</E><v>0.3</v>", ""),
    "C_iso_truss": ("isotropic elastic",
                     "      <density>1</density><E>260</E><v>0.3</v>", ""),
}

for name, (mt, params, extra) in VARIANTS.items():
    f = WORK / f"{name}.feb"
    f.write_text(HEAD.format(mt=mt, params=params, extra=extra), encoding="latin-1")
    r = subprocess.run([str(FEBIO), "-i", str(f)], cwd=str(WORK),
                       capture_output=True, text=True, timeout=300)
    o = (r.stdout or "") + (r.stderr or "")
    ok = "NORMAL TERMINATION" in o
    err = [l.strip() for l in o.split("\n")
           if "ERROR" in l.upper() or "invalid" in l.lower()
           or "unrecognized" in l.lower() or "must" in l.lower()][:4]
    print(f"\n--- {name} (type={mt}) rc={r.returncode} "
          f"{'NORMAL TERMINATION' if ok else 'FAILED'} ---")
    for e in err:
        print("   ", e)
    for l in o.split("\n"):
        if "nlx" in l.lower() or "strain" in l.lower() or "displacement" in l.lower():
            print("   >", l.strip()[:120])
