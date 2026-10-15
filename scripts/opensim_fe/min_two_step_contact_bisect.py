# -*- coding: utf-8 -*-
"""bisection 第二轮：**加 `<Contact>`**（+ 第二块 + SurfacePair）。

第二轮动机：`min_two_step_bisect.py` 第一轮已排除 `<LoadData>`/`lc`/三向 prescribe
（力步初残差仍 38.9 / 657 ⇒ 载荷有效）。跟骨 deck 剩下的重量级差异就是 **`<Contact>`**
（单边约束 + `pad_floor` 第三域 + `<SurfacePair>`）。

模型：两块 hex8，A 在 x∈[-2,2]、B 在 x∈[2.1,3.1]（**初始间隙 0.1**），
A 右面 ⇄ B 左面 组成 `<SurfacePair>`；`<Contact>` 放**全局**（§6.1.1：各步都生效）。
Step 1 用 prescribed displacement 把 A 压过去（+0.15 > 0.1 间隙）；Step 2 换力载荷。

**判据**：力步初残差是否仍为载荷量级（~1e1）。
用法::

    .venv\\Scripts\\python.exe scripts\\opensim_fe\\min_two_step_contact_bisect.py
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

OUT = ROOT / "temp" / "opensim_fe" / "min2step_contact"

TMPL = """<?xml version="1.0" encoding="ISO-8859-1"?>
<febio_spec version="4.0">
  <Module type="solid"/>
  <Material>
    <material id="1" name="MA" type="neo-Hookean">
      <density>1.0</density><E>1000</E><v>0.3</v>
    </material>
    <material id="2" name="MB" type="neo-Hookean">
      <density>1.0</density><E>100000</E><v>0.3</v>
    </material>
  </Material>
  <Mesh>
    <Nodes>
      <node id="1">-2.0,-0.5, 0.0</node>
      <node id="2">-2.0,-0.5, 1.0</node>
      <node id="3">-2.0, 0.5, 0.0</node>
      <node id="4">-2.0, 0.5, 1.0</node>
      <node id="5"> 2.0,-0.5, 0.0</node>
      <node id="6"> 2.0,-0.5, 1.0</node>
      <node id="7"> 2.0, 0.5, 0.0</node>
      <node id="8"> 2.0, 0.5, 1.0</node>
      <node id="9"> 2.1,-0.5, 0.0</node>
      <node id="10"> 2.1,-0.5, 1.0</node>
      <node id="11"> 2.1, 0.5, 0.0</node>
      <node id="12"> 2.1, 0.5, 1.0</node>
      <node id="13"> 3.1,-0.5, 0.0</node>
      <node id="14"> 3.1,-0.5, 1.0</node>
      <node id="15"> 3.1, 0.5, 0.0</node>
      <node id="16"> 3.1, 0.5, 1.0</node>
    </Nodes>
    <Elements type="hex8" name="PartA"><elem id="1">1,5,7,3,2,6,8,4</elem></Elements>
    <Elements type="hex8" name="PartB"><elem id="2">9,13,15,11,10,14,16,12</elem></Elements>
    <NodeSet name="set1">1,2,3,4</NodeSet>
    <NodeSet name="setB">13,14,15,16</NodeSet>
    <NodeSet name="setAright">5,6,7,8</NodeSet>
    <Surface name="Aright"><quad4 id="1">5,6,8,7</quad4></Surface>
    <Surface name="Bleft"><quad4 id="1">{BLEFT}</quad4></Surface>
    <SurfacePair name="pairAB">
      <primary>{PRIM}</primary><secondary>{SEC}</secondary>
    </SurfacePair>
  </Mesh>
  <MeshDomains>
    <SolidDomain name="PartA" mat="MA"/>
    <SolidDomain name="PartB" mat="MB"/>
  </MeshDomains>
  <Boundary>
    <bc type="zero displacement" node_set="set1">
      <x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof>
    </bc>
  </Boundary>
  <Contact>
    <contact name="pad" surface_pair="pairAB" type="sliding-elastic">
      <laugon>PENALTY</laugon>
      <penalty>0.1</penalty>
      <auto_penalty>0</auto_penalty>
      <tolerance>0.005</tolerance>
      <search_radius>1.0</search_radius>
      <fric_coeff>0.6</fric_coeff>
    </contact>
  </Contact>
  <Step>
    <step id="1" name="pretouch">
      <Control><time_steps>2</time_steps><step_size>0.5</step_size><solver/></Control>
      <Boundary>
        <bc name="push" type="prescribed displacement" node_set="setAright">
          <dof>x</dof><value>0.15</value>
        </bc>
      </Boundary>
    </step>
    <step id="2" name="force">
      <Control><time_steps>2</time_steps><step_size>0.5</step_size><solver/></Control>
      <Loads>
        <surface_load name="p" type="pressure" surface="Aright">
          <pressure>1.0</pressure>
        </surface_load>
      </Loads>
    </step>
  </Step>
  <Output><plotfile type="febio" file="{TAG}.xplt"><var type="displacement"/></plotfile></Output>
</febio_spec>
"""

VARIANTS = {
    # B 左面 (x=2.1) 的 quad4 节点序；primary/secondary 两种配对
    "c0_A_prim":  dict(BLEFT="9,10,12,11", PRIM="Aright", SEC="Bleft"),
    "c1_B_prim":  dict(BLEFT="9,10,12,11", PRIM="Bleft", SEC="Aright"),
    "c2_rev":     dict(BLEFT="9,11,12,10", PRIM="Aright", SEC="Bleft"),
}


def main() -> int:
    from climbing.coupling.febio_run import find_febio

    OUT.mkdir(parents=True, exist_ok=True)
    exe = find_febio()
    print('%-12s %-4s %-18s %-14s %s' % ('variant', 'rc', '力步初残差', 'negJ', '结论'))
    print('-' * 82)
    for tag, spec in VARIANTS.items():
        feb = OUT / ('%s.feb' % tag)
        for s in ('.xplt', '.log'):
            p = feb.with_suffix(s)
            if p.exists():
                p.unlink()
        feb.write_bytes(TMPL.format(TAG=tag, **spec).encode("iso-8859-1"))
        p = subprocess.run([str(exe), "-i", str(feb)], capture_output=True,
                           text=True, cwd=str(OUT), timeout=300)
        log = feb.with_suffix(".log")
        t = log.read_text(encoding="utf-8", errors="replace") if log.is_file() else ""
        blocks = list(re.finditer(r'beginning time step 1\s*:\s*([\d.]+)', t))
        s2 = '(no 2nd block)'
        if len(blocks) >= 2:
            r = re.search(r'residual\s+([\d.eE+-]+)\s+([\d.eE+-]+)',
                          t[blocks[1].end():])
            if r:
                s2 = r.group(1)
        nj = re.findall(r'(\d+)\s*negative\s+jacobians?', t, re.I)
        try:
            mag = float(s2) > 1e-6
        except (TypeError, ValueError):
            mag = False
        ok = p.returncode == 0
        print('%-12s %-4d %-18s %-14s %s' % (
            tag, p.returncode, s2, nj[-1] if nj else '-',
            ('✅ 载荷有效' if (ok and mag) else
             ('⚠️ 通过但载荷塌 0' if ok else '❌ 失败'))))
    print()
    print('判据：力步初残差应 ~1e1。塌成 ~1e-16 ⇒ **Contact 就是根因**。')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
