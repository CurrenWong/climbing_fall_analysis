"""最小验证 FEBio 刚体机制：rigid-solid 域 + rigid BC + rigid spring + rigid_rotation。

设计：
  body A = 1 个 rigid hex8 立方体（要移动的"距骨"）
  body G = 1 个节点（"地面"），全约束
  A、G 之间一根 rigid spring，k=500 N/mm，沿 z
  A 上施加 rigid_force 100 N 沿 z
判据： A 的 z 位移 = F/k = 100/500 = 0.2 mm
      （若 <rigid> 机制不通 ⇒ 报错或位移不是 0.2）
"""
import re
import subprocess
from pathlib import Path

W = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo")
FEBIO = r"D:\Program\FEBioStudio\bin\febio4.exe"

TMPL = """<?xml version="1.0" encoding="ISO-8859-1"?>
<febio_spec version="4.0">
  <Module type="solid"/>
  <Control>
    <analysis>STATIC</analysis><time_steps>1</time_steps><step_size>1.000000</step_size>
    <solver type="solid">
      <linear_solver type="pardiso"/>
      <max_refs>50</max_refs>
      <diverge_reform>1</diverge_reform><reform_each_time_step>1</reform_each_time_step>
      <dtol>0.01</dtol><etol>0.001</etol><rtol>0</rtol><lstol>0.9</lstol>
      <min_residual>1e-20</min_residual><rhoi>-2</rhoi>
      <qn_method type="BFGS"><max_ups>10</max_ups></qn_method>
    </solver>
    <time_stepper>
      <dtmin>0.0001</dtmin><dtmax>0.05</dtmax>
      <max_retries>10</max_retries><opt_iter>10</opt_iter>
    </time_stepper>
  </Control>
  <Material>
    <material id="1" type="rigid body" name="bodyA">
      <density>1e-9</density>
      <E>12000</E><v>0.3</v>
    </material>
    <material id="2" type="rigid body" name="bodyG">
      <density>1e-9</density>
      <E>12000</E><v>0.3</v>
    </material>
  </Material>
  <Mesh>
    <Nodes>
      <node id="1">0,0,0</node>
      <node id="2">1,0,0</node>
      <node id="3">1,1,0</node>
      <node id="4">0,1,0</node>
      <node id="5">0,0,1</node>
      <node id="6">1,0,1</node>
      <node id="7">1,1,1</node>
      <node id="8">0,1,1</node>
      <node id="9">0,0,-1</node>
    </Nodes>
    <Elements type="hex8" name="partA">
      <elem id="1">1,2,3,4,5,6,7,8</elem>
    </Elements>
    <NodeSet name="nsA">1,2,3,4,5,6,7,8</NodeSet>
    <NodeSet name="nsG">9</NodeSet>
  </Mesh>
  <MeshDomains>
    <SolidDomain name="partA" type="rigid-solid" mat="bodyA"/>
  </MeshDomains>
  <LoadData>
    <load_controller id="1" type="loadcurve">
      <interpolate>LINEAR</interpolate>
      <points>
        <point>0,0</point>
        <point>1,1</point>
      </points>
    </load_controller>
  </LoadData>
  <Boundary>
    <bc type="rigid" node_set="nsA"><rb>bodyA</rb></bc>
    <bc type="rigid" node_set="nsG"><rb>bodyG</rb></bc>
  </Boundary>
  <Rigid>
    <rigid_bc type="rigid_fixed">
      <rb>bodyG</rb>
      <Rx_dof>1</Rx_dof><Ry_dof>1</Ry_dof><Rz_dof>1</Rz_dof>
      <Ru_dof>1</Ru_dof><Rv_dof>1</Rv_dof><Rw_dof>1</Rw_dof>
    </rigid_bc>
    <rigid_bc type="rigid_fixed">
      <rb>bodyA</rb>
      <Ru_dof>{FIXRU}</Ru_dof><Rv_dof>1</Rv_dof><Rw_dof>{FIXRW}</Rw_dof>
      <Rx_dof>1</Rx_dof><Ry_dof>1</Ry_dof>
    </rigid_bc>
    <rigid_connector type="rigid spring">
      <body_a>1</body_a>
      <body_b>2</body_b>
      <insertion_a>0.5,0.5,0.5</insertion_a>
      <insertion_b>0,0,-1</insertion_b>
      <k>500</k>
    </rigid_connector>
    <rigid_load type="rigid_force">
      <rb>bodyA</rb>
      <dof>Rz</dof>
      <value lc="1">100</value>
    </rigid_load>
  </Rigid>
  <Output>
    <plotfile type="febio"><var type="displacement"/>
      <var type="reaction forces"/></plotfile>
    <logfile file="log_rigid.txt">
      <node_data data="ux;uy;uz" name="dispA" file="disp_rigid.txt" delim=",">1,5</node_data>
    </logfile>
  </Output>
</febio_spec>
"""

for tag, fixru, fixrw in (("locked_rot", "1", "1"), ("free_rot", "0", "0")):
    f = W / f"rigid_min_{tag}.feb"
    f.write_text(TMPL.replace("{FIXRU}", fixru).replace("{FIXRW}", fixrw),
                 encoding="utf-8")
    r = subprocess.run([FEBIO, "-i", f.name, "-o", f"{f.stem}.xplt"],
                       capture_output=True, text=True, cwd=str(W), timeout=300)
    out = (r.stdout or "") + (r.stderr or "")
    ok = ("N O R M A L" in out and "T E R M I N A T I O N" in out) or \
         "NORMAL TERMINATION" in out
    print(f"\n{'='*70}\n### {tag}  (Ru fixed={fixru})\n{'='*70}")
    print(f"  收敛: {'YES ✅' if ok else 'NO ❌'}")
    d = W / "disp_rigid.txt"
    if d.exists():
        rows = [l for l in d.read_text(encoding="utf-8", errors="replace").splitlines()
                if l.strip() and not l.startswith("*")]
        print(f"  最后一组位移读数（期望 uz = 100/500 = 0.2 mm）:")
        for l in rows[-2:]:
            print(f"    {l.strip()}")
    if not ok:
        for ln in out.splitlines():
            if any(x in ln.lower() for x in ("error", "tag ", "invalid",
                                             "unrecognized", "failed")):
                print("   " + ln.strip()[:110])
