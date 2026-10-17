"""最小验证 FEBio <Discrete> 线性弹簧是否被真正加载。

设计：立方体底面全固定，顶面 4 节点 z 方向指定位移 δ=-0.1。
      额外一根弹簧从顶面节点 5 连到固定的哑节点 9（沿 z）。
判据：**顶面 z 支反力之差** 必须 ≈ (k2-k1)*δ。
      若 <Discrete> 被忽略 ⇒ 差值 = 0 ⇒ 语法无效。
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
      <diverge_reform>1</diverge_reform>
      <reform_each_time_step>1</reform_each_time_step>
      <dtol>0.01</dtol>
      <etol>0.001</etol>
      <rtol>0</rtol>
      <lstol>0.9</lstol>
      <min_residual>1e-20</min_residual>
      <rhoi>-2</rhoi>
      <qn_method type="BFGS"><max_ups>10</max_ups></qn_method>
    </solver>
    <time_stepper>
      <dtmin>0.0001</dtmin><dtmax>0.05</dtmax>
      <max_retries>10</max_retries><opt_iter>10</opt_iter>
    </time_stepper>
  </Control>
  <Material>
    <material id="1" name="m1" type="isotropic elastic">
      <density>1e-9</density>
      <E>1000</E>
      <v>0.3</v>
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
      <node id="9">0,0,2</node>
    </Nodes>
    <Elements type="hex8" name="cube">
      <elem id="1">1,2,3,4,5,6,7,8</elem>
    </Elements>
    <NodeSet name="bottom">1,2,3,4</NodeSet>
    <NodeSet name="top">5,6,7,8</NodeSet>
    <NodeSet name="dummy">9</NodeSet>
    <DiscreteSet name="springs">
      <delem>5,9</delem>
    </DiscreteSet>
  </Mesh>
  <MeshDomains>
    <SolidDomain name="cube" mat="m1"/>
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
  <Discrete>
    <discrete_material id="1" type="linear spring">
      <E>{K}</E>
    </discrete_material>
    <discrete dmat="1" discrete_set="springs"/>
  </Discrete>
  <Boundary>
    <bc name="fixbot" type="zero displacement" node_set="bottom">
      <x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof>
    </bc>
    <bc name="fixdum" type="zero displacement" node_set="dummy">
      <x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof>
    </bc>
    <bc name="pull" type="prescribed displacement" node_set="top">
      <dof>z</dof><value lc="1">-0.1</value><relative>0</relative>
    </bc>
  </Boundary>
  <Output>
    <plotfile type="febio"><var type="displacement"/><var type="reaction forces"/></plotfile>
    <logfile file="log_springtest.txt">
      <node_data data="Rz" name="reac" file="reac_springtest.txt" delim=",">5,6,7,8</node_data>
    </logfile>
  </Output>
</febio_spec>
"""

print(f"{'k':>8} {'rc':>4} {'Σ顶面 Rz (N)':>14}  {'log':>6}")
print("-" * 44)
res = {}
for k in (0.001, 1.0, 100.0, 200.0):
    f = W / f"spring_min_{k}.feb"
    f.write_text(TMPL.replace("{K}", str(k)), encoding="utf-8")
    for junk in ("reac_springtest.txt",):
        p = W / junk
        if p.exists():
            p.unlink()
    r = subprocess.run([FEBIO, "-i", f.name, "-o", f"{f.stem}.xplt"],
                       capture_output=True, text=True, cwd=str(W), timeout=300)
    out = (r.stdout or "") + (r.stderr or "")
    rcf = re.search(r"FEBio rc=(-?\d+)", out)
    norm = "NORMAL" in out and "TERMINATION" in out
    reac = None
    rp = W / "reac_springtest.txt"
    if rp.exists():
        vals = []
        for ln in rp.read_text(encoding="utf-8", errors="replace").splitlines():
            for x in re.findall(r"[-\d.eE+]+", ln):
                try:
                    vals.append(float(x))
                except ValueError:
                    pass
        reac = sum(vals) if vals else None
    res[k] = reac
    print(f"{k:>8} {('OK' if norm else 'ERR'):>4} "
          f"{(f'{reac:.4f}' if reac is not None else '-'):>14}  "
          f"{'有' if rp.exists() else '无'}")

d = res.get(100.0, 0) - res.get(1.0, 0)
print(f"\nΔF(100 vs 1) = {d:.4f} N    期望 = Δk*δ = 99*0.1 = 9.9 N")
print("⇒", "**弹簧生效 ✅**" if abs(d - 9.9) < 0.5 else "⚠ 不符，语法可能未被采纳")
