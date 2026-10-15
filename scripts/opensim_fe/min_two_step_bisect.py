# -*- coding: utf-8 -*-
"""bisection：把跟骨 2 步 deck 的特征逐个加进最小复现器，找「加到哪一步开始坏」。

`min_two_step_load_test.py` 的 v0（`<Loads>` 在 step 内）**通过**（Step 2 初残差 38.9）。
本件在其上逐项加跟骨 deck 才有的特征，每项单变量：

| 变体 | 加什么 | 对应跟骨 deck 的哪一处 |
|---|---|---|
| `v0_base` | — | 最小基线（已知 ✅） |
| `v1_ldata` | 顶层 `<LoadData>`（lc=1 曲线 0→1） | 跟骨 deck 有（为 prescribe 的 `lc="1"`） |
| `v2_3dof` | Step 1 prescribe **三向** | 跟骨 deck 三向全控 |
| `v3_lc_load` | 压力载荷挂 `lc="1"` | （对照：跟骨 deck 的载荷**没有** lc） |
| `v4_ldata_3dof` | v1 + v2 | 跟骨 deck 的真实组合 |

**判据**：Step 2 的「初始残差」是否仍为**载荷量级**（~1e1），还是塌成 ~1e-16。
用法::

    .venv\\Scripts\\python.exe scripts\\opensim_fe\\min_two_step_bisect.py
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

OUT = ROOT / "temp" / "opensim_fe" / "min2step_bisect"

BODY = """<?xml version="1.0" encoding="ISO-8859-1"?>
<febio_spec version="4.0">
  <Module type="solid"/>
  <Material>
    <material id="1" name="M" type="neo-Hookean">
      <density>1.0</density><E>1000</E><v>0.3</v>
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
    <Elements type="hex8" name="Part1"><elem id="1">1,5,7,3,2,6,8,4</elem></Elements>
    <NodeSet name="set1">1,2,3,4</NodeSet>
    <NodeSet name="set4">5,6,7,8</NodeSet>
    <Surface name="right"><quad4 id="1">5,6,8,7</quad4></Surface>
  </Mesh>
  <MeshDomains><SolidDomain name="Part1" mat="M"/></MeshDomains>
  <Boundary>
    <bc type="zero displacement" node_set="set1">
      <x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof>
    </bc>
  </Boundary>
{LDATA}  <Step>
    <step id="1" name="prescribe">
      <Control><time_steps>2</time_steps><step_size>0.5</step_size><solver/></Control>
      <Boundary>
{PRESCRIBE}      </Boundary>
    </step>
    <step id="2" name="force">
      <Control><time_steps>2</time_steps><step_size>0.5</step_size><solver/></Control>
      <Loads>
        <surface_load name="p" type="pressure" surface="right"{LC}>
          <pressure>1.0</pressure>
        </surface_load>
      </Loads>
    </step>
  </Step>
  <Output><plotfile type="febio" file="{TAG}.xplt"><var type="displacement"/></plotfile></Output>
</febio_spec>
"""

LDATA = """  <LoadData>
    <load_controller id="1" type="loadcurve">
      <interpolate>LINEAR</interpolate><extend>CONSTANT</extend>
      <points><point>0,0</point><point>1,1</point></points>
    </load_controller>
  </LoadData>
"""

BC1 = ('        <bc name="push" type="prescribed displacement" node_set="set4">\n'
       '          <dof>{d}</dof><value>0.05</value>\n'
       '        </bc>\n')

VARIANTS = {
    "v0_base":        dict(ldata=False, dofs="x", lc=False),
    "v1_ldata":       dict(ldata=True,  dofs="x", lc=False),
    "v2_3dof":        dict(ldata=False, dofs="xyz", lc=False),
    "v3_lc_load":     dict(ldata=True,  dofs="x", lc=True),
    "v4_ldata_3dof":  dict(ldata=True,  dofs="xyz", lc=False),
}


def _deck(tag: str, spec: dict) -> str:
    pre = "".join(BC1.format(d=d) for d in spec["dofs"])
    return BODY.format(TAG=tag,
                       LDATA=(LDATA if spec["ldata"] else ""),
                       PRESCRIBE=pre,
                       LC=(' lc="1"' if spec["lc"] else ""))


def main() -> int:
    from climbing.coupling.febio_run import find_febio

    OUT.mkdir(parents=True, exist_ok=True)
    exe = find_febio()
    print('%-16s %-4s %-22s %s' % ('variant', 'rc', 'Step2 初残差', '结论'))
    print('-' * 78)
    for tag, spec in VARIANTS.items():
        feb = OUT / ('%s.feb' % tag)
        for s in ('.xplt', '.log'):
            p = feb.with_suffix(s)
            if p.exists():
                p.unlink()
        feb.write_bytes(_deck(tag, spec).encode("iso-8859-1"))
        p = subprocess.run([str(exe), "-i", str(feb)], capture_output=True,
                           text=True, cwd=str(OUT), timeout=300)
        log = feb.with_suffix(".log")
        t = log.read_text(encoding="utf-8", errors="replace") if log.is_file() else ""
        # ⚠️ 日志里 **`beginning time step 1` 会出现两次** —— FEBio 对**每个 `<step>`
        # 各自从 `time step 1` 开始编号**（累加时间照走：1.5 → 2）。故：
        #   * 第 1 次 `time step 1` 块 ⇒ prescribe 步
        #   * 第 2 次 `time step 1` 块 ⇒ **force 步**（取它的第一条 residual 第一列）
        blocks = list(re.finditer(r'beginning time step 1\s*:\s*([\d.]+)', t))
        s2_init = '(no 2nd step-1 block)'
        if len(blocks) >= 2:
            r = re.search(r'residual\s+([\d.eE+-]+)\s+([\d.eE+-]+)',
                          t[blocks[1].end():])
            if r:
                s2_init = r.group(1)
        ok = p.returncode == 0
        try:
            mag = float(s2_init) > 1e-6
        except (TypeError, ValueError):
            mag = False
        print('%-16s %-4d %-22s %s' % (
            tag, p.returncode, s2_init,
            ('✅ 通过，载荷有效' if (ok and mag)
             else ('⚠️ 通过但**载荷塌为 0**' if ok else '❌ 失败'))))
    print()
    print('判据：Step2 初残差应 ~1e1（载荷量级）。塌成 ~1e-16 ⇒ 该增量就是根因。')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
