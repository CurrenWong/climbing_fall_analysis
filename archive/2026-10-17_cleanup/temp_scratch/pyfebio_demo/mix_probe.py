"""最小实验：杆单元(BeamDomain) 与实体单元共用节点 —— FEBio 是否支持？"""
import subprocess
from pathlib import Path

W = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\mix_probe")
W.mkdir(parents=True, exist_ok=True)
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")

MM = '''    <material id="1" name="solid" type="isotropic elastic">
      <density>2e-9</density><E>7300</E><v>0.3</v>
    </material>
    <material id="2" name="lig" type="isotropic elastic">
      <density>1.1e-9</density><E>5200</E><v>0.4</v>
    </material>'''

HEAD = '''<?xml version="1.0" encoding="ISO-8859-1"?>
<febio_spec version="4.0">
  <Module type="solid"/>
  <Globals><Constants><T>0</T><R>0</R><Fc>0</Fc></Constants></Globals>
  <Material>
{mat}
  </Material>
  <Mesh>
    <Nodes name="all">
{elem_nodes}
{extra_nodes}
    </Nodes>
    <Elements type="hex8" name="blk">
      <elem id="1">1,2,3,4,5,6,7,8</elem>
    </Elements>
{elem_lig}
    <NodeSet name="fix">1,2,3,4</NodeSet>
    <NodeSet name="pull">{pull}</NodeSet>
  </Mesh>
  <MeshDomains>
    <SolidDomain name="blk" mat="solid"/>
{dom_lig}
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
      <dof>x</dof><value lc="1">-0.5</value><relative>0</relative>
    </bc>
  </Boundary>
  <Output>
    <plotfile type="febio"><var type="displacement"/></plotfile>
  </Output>
  <Control>
    <analysis>STATIC</analysis><time_steps>2</time_steps><step_size>0.5</step_size>
    <solver type="solid"><linear_solver type="pardiso"/></solver>
  </Control>
</febio_spec>
'''

# 一个 10x10x10 的六面体：节点 1-8
HEX = "\n".join(f"      <node id=\"{i}\">{p}</node>" for i, p in enumerate([
    (0, 0, 0), (10, 0, 0), (10, 10, 0), (0, 10, 0),
    (0, 0, 10), (10, 0, 10), (10, 10, 10), (0, 10, 10)], 1))

CASES = {
    # A: 杆两端都接在实体节点上（共享节点）
    "A_shared_both": dict(
        extra_nodes="", elem_lig='''    <Elements type="line2" name="lt">
      <elem id="2">2,4</elem>
    </Elements>''',
        dom_lig='    <BeamDomain name="lt" mat="lig" type="elastic-truss"/>',
        pull="4"),
    # B: 杆有一端是独立节点
    "B_one_free": dict(
        extra_nodes='      <node id="9">20,5,5</node>',
        elem_lig='''    <Elements type="line2" name="lt">
      <elem id="2">2,9</elem>
    </Elements>''',
        dom_lig='    <BeamDomain name="lt" mat="lig" type="elastic-truss"/>',
        pull="9"),
    # C: 无杆，只实体（对照）
    "C_solid_only": dict(
        extra_nodes="", elem_lig="",
        dom_lig="", pull="7"),
}

for name, kw in CASES.items():
    t = HEAD.format(mat=MM, elem_nodes=HEX, **kw)
    f = W / f"{name}.feb"
    f.write_text(t, encoding="latin-1")
    r = subprocess.run([str(FEBIO), "-i", str(f)], cwd=str(W),
                       capture_output=True, text=True, timeout=300)
    o = (r.stdout or "") + (r.stderr or "")
    ok = "NORMAL" in o.replace(" ", "")
    print(f"\n=== {name} rc={r.returncode} {'OK' if ok else 'FAIL'} ===")
    lines = [l.strip() for l in o.split("\n") if l.strip()]
    for l in lines:
        if "ERROR" in l.upper() or "invalid" in l.lower() \
           or "unrecognized" in l.lower() or "SUCCESS" in l:
            print("   ", l[:130])
    if not any("ERROR" in l.upper() for l in lines):
        print("   尾:", lines[-1][:100] if lines else "")
