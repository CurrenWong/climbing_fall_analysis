"""二分：truss 节点接在 tet4 上 / 多材料 / Surface / SurfacePair 是否触发崩溃。"""
import subprocess
from pathlib import Path

W = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\mix_probe")
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")

M2 = '''    <material id="1" name="solid" type="isotropic elastic">
      <density>2e-9</density><E>7300</E><v>0.3</v>
    </material>
    <material id="2" name="lig" type="isotropic elastic">
      <density>1.1e-9</density><E>5200</E><v>0.4</v>
    </material>
    <material id="3" name="cart" type="isotropic elastic">
      <density>1.1e-9</density><E>260</E><v>0.4</v>
    </material>'''


def build(etype, extra="", extra_dom="", extra_mat="", lig_elem="",
          lig_dom="", nodes_extra=""):
    if etype == "hex8":
        nodes = "".join(
            f'      <node id="{i}">{p[0]},{p[1]},{p[2]}</node>'
            for i, p in enumerate([(0, 0, 0), (10, 0, 0), (10, 10, 0),
                                   (0, 10, 0), (0, 0, 10), (10, 0, 10),
                                   (10, 10, 10), (0, 10, 10)], 1))
        elems = '<elem id="1">1,2,3,4,5,6,7,8</elem>'
        fix = "1,2,3,4"
    else:
        nodes = "".join(
            f'      <node id="{i}">{p[0]},{p[1]},{p[2]}</node>'
            for i, p in enumerate([(0, 0, 0), (10, 0, 0), (0, 10, 0),
                                   (0, 0, 10)], 1))
        elems = '<elem id="1">1,2,3,4</elem>'
        fix = "1,2,3"
    return f'''<?xml version="1.0" encoding="ISO-8859-1"?>
<febio_spec version="4.0">
  <Module type="solid"/>
  <Globals><Constants><T>0</T><R>0</R><Fc>0</Fc></Constants></Globals>
  <Material>
{M2}
{extra_mat}
  </Material>
  <Mesh>
    <Nodes name="all">
{nodes}
{nodes_extra}
    </Nodes>
    <Elements type="{etype}" name="blk">
      {elems}
    </Elements>
{lig_elem}
    <NodeSet name="fix">{fix}</NodeSet>
    <NodeSet name="pull">2</NodeSet>
    <Surface name="s_top">
      <tri3 id="1">1,2,3</tri3>
    </Surface>
    <SurfacePair name="sp">
      <primary>s_top</primary>
      <secondary>s_top</secondary>
    </SurfacePair>
    {extra}
  </Mesh>
  <MeshDomains>
    <SolidDomain name="blk" mat="solid"/>
{extra_dom}
{lig_dom}
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
  </Boundary>
  <Contact>
    <contact name="c" surface_pair="sp" type="sliding-elastic">
      <laugon>PENALTY</laugon><penalty>1</penalty><auto_penalty>1</auto_penalty>
    </contact>
  </Contact>
  <Output>
    <plotfile type="febio"><var type="displacement"/></plotfile>
  </Output>
  <Control>
    <analysis>STATIC</analysis><time_steps>2</time_steps><step_size>0.5</step_size>
    <solver type="solid"><linear_solver type="pardiso"/></solver>
  </Control>
</febio_spec>
'''


LIGE_TET = '''    <Elements type="line2" name="lt">
      <elem id="9">1,4</elem>
    </Elements>'''
LIGD_TET = '    <BeamDomain name="lt" mat="lig" type="elastic-truss"/>'

CASES = {
    "T1_tet_plus_truss": build("tet4", lig_elem=LIGE_TET, lig_dom=LIGD_TET),
    "T2_tet_only": build("tet4"),
    "T3_tet_truss_nocontact": build("tet4", lig_elem=LIGE_TET, lig_dom=LIGD_TET),
}

for name, t in CASES.items():
    if name == "T3_tet_truss_nocontact":
        i0 = t.index("  <Contact>")
        i1 = t.index("</Contact>") + len("</Contact>") + 1
        t = t[:i0] + t[i1:]
    f = W / f"{name}.feb"
    f.write_text(t, encoding="latin-1")
    r = subprocess.run([str(FEBIO), "-i", str(f)], cwd=str(W),
                       capture_output=True, text=True, timeout=300)
    o = (r.stdout or "") + (r.stderr or "")
    ok = "NORMAL" in o.replace(" ", "")
    print(f"\n=== {name} rc={r.returncode} {'OK' if ok else 'FAIL'} ===")
    for l in [x.strip() for x in o.split("\n") if x.strip()]:
        if "ERROR" in l.upper() or "invalid" in l.lower() \
           or "unrecognized" in l.lower():
            print("   ", l[:130])
