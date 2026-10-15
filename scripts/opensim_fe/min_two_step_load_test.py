# -*- coding: utf-8 -*-
"""最小复现器：**`<Loads>` 放进 `<step>` 到底生不生效？**（手册 §6.2 单单元模型）

诊断背景（`two_step_pretouch_probe.py`）：跟骨 2 步 deck 里 Step 1（prescribed displacement）
成功、**Step 2（力载荷）残差恒 ~1e-16** ⇒ 疑「载荷未进求解」。在跟骨模型上换 4 种放置都无效
（step 内 / 全局 / 删顶层 Control / 加 `<solver><min_residual>`）⇒ 需**最小体系**隔离。

本件只用**手册 §6.2 的 8 节点单单元模型**，把它扩成两步：
  Step 1：prescribed displacement（x）+  **pressure 载荷（已知正对照，小量）**
  Step 2：**只** 放一个 pressure 载荷

判据（**已知正对照**，纪律 4）：
  * Step 2 的 `residual INITIAL` 必须 ≈ 载荷量级（~1），**不是** ~1e-16
  * 若 Step 2 的 residual 也是 ~1e-16 ⇒ **`<Loads>` 在 step 内确实不生效** ⇒ 结构问题
用法::

    .venv\\Scripts\\python.exe scripts\\opensim_fe\\min_two_step_load_test.py
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

OUT = ROOT / "temp" / "opensim_fe" / "min2step"

HEAD = """<?xml version="1.0" encoding="ISO-8859-1"?>
<febio_spec version="4.0">
  <Module type="solid"/>
  <Material>
    <material id="1" name="M" type="neo-Hookean">
      <density>1.0</density>
      <E>1000</E>
      <v>0.3</v>
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
    </Nodes>
    <Elements type="hex8" name="Part1">
      <elem id="1">1,5,7,3,2,6,8,4</elem>
    </Elements>
    <NodeSet name="set1">1,2,3,4</NodeSet>
    <NodeSet name="set4">5,6,7,8</NodeSet>
    <Surface name="right">
      <quad4 id="1">5,6,8,7</quad4>
    </Surface>
  </Mesh>
  <MeshDomains>
    <SolidDomain name="Part1" mat="M"/>
  </MeshDomains>
  <Boundary>
    <bc type="zero displacement" node_set="set1">
      <x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof>
    </bc>
  </Boundary>
  <Step>
    <step id="1" name="prescribe">
      <Control><time_steps>2</time_steps><step_size>0.5</step_size>
        <solver/></Control>
      <Boundary>
        <bc name="push" type="prescribed displacement" node_set="set4">
          <dof>x</dof><value>0.05</value>
        </bc>
      </Boundary>
    </step>
    <step id="2" name="force">
      <Control><time_steps>2</time_steps><step_size>0.5</step_size>
        <solver/></Control>
      <Loads>
        <surface_load name="p" type="pressure" surface="right">
          <pressure>1.0</pressure>
        </surface_load>
      </Loads>
    </step>
  </Step>
  <Output>
    <plotfile type="febio" file="%s.xplt"><var type="displacement"/></plotfile>
  </Output>
</febio_spec>
"""


def main() -> int:
    from climbing.coupling.febio_run import find_febio

    OUT.mkdir(parents=True, exist_ok=True)
    feb = OUT / "min2step.feb"
    for s in (".xplt", ".log"):
        p = feb.with_suffix(s)
        if p.exists():
            p.unlink()
    feb.write_bytes((HEAD % "min2step").encode("iso-8859-1"))
    d = feb.read_bytes().decode("iso-8859-1", errors="replace")
    print('[deck]', feb)
    print('  [闸门3] steps=%d  prescribed=%d  Loads=%d  surfaces=%d' % (
        len(re.findall(r'<step\b', d)),
        len(re.findall(r'prescribed displacement', d)),
        len(re.findall(r'<Loads\b', d)),
        len(re.findall(r'<Surface\b', d))))

    exe = find_febio()
    proc = subprocess.run([str(exe), "-i", str(feb)], capture_output=True,
                          text=True, cwd=str(OUT), timeout=300)
    print('  [solve] rc=%d' % proc.returncode)
    log = feb.with_suffix(".log")
    t = log.read_text(encoding="utf-8", errors="replace") if log.is_file() else ""
    if not t:
        print((proc.stdout or "")[-1800:])
        print((proc.stderr or "")[-800:])
        return 1
    print('  --- 关键行 ---')
    for l in t.split('\n'):
        if re.search(r'(beginning time step|residual|failed to converge|'
                     r'negative|Number of time steps completed|NORMAL TERMINATION|'
                     r'ERROR)', l, re.I):
            s = l.strip()
            if s:
                print('    ', s[:120])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
