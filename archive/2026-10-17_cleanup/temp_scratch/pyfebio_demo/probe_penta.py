"""最小实验：确定 FEBio 4.13 的 penta6 节点序约定。

底面三角形 (0,0,0) (1,0,0) (0,1,0)，顶面 (0,0,1) (1,0,1) (0,1,1)。
分别按若干排列写成单单元 .feb，跑 FEBio 看是否报
"Negative jacobian detected during mesh initialization"。
"""
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
TMP = ROOT / "temp/pyfebio_demo/penta_probe"
TMP.mkdir(parents=True, exist_ok=True)

COORD = {
    1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0), 3: (0.0, 1.0, 0.0),
    4: (0.0, 0.0, 1.0), 5: (1.0, 0.0, 1.0), 6: (0.0, 1.0, 1.0),
}

VARIANTS = {
    "A_0,1,2,3,4,5 (底+顶同序)": [1, 2, 3, 4, 5, 6],
    "B_0,1,2,6,5,4 (顶反序)":    [1, 2, 3, 6, 5, 4],
    "C_0,1,2,5,6,4":             [1, 2, 3, 5, 6, 4],
    "D_0,1,2,4,6,5":             [1, 2, 3, 4, 6, 5],
    "E_0,1,2,4,5,6 (底+顶同)":   [1, 2, 3, 4, 5, 6],
    "F_1,0,2,4,5,6":             [2, 1, 3, 4, 5, 6],
}


def build(conn):
    nodes = "\n".join(f'      <node id="{i}">{x},{y},{z}</node>'
                      for i, (x, y, z) in COORD.items())
    return f"""<?xml version="1.0" encoding="utf-8"?>
<febio_spec version="4.0">
  <Module type="solid"/>
  <Material>
    <material id="1" name="m" type="isotropic elastic">
      <density>1e-9</density>
      <E>10.0</E>
      <v>0.4</v>
    </material>
  </Material>
  <Mesh>
    <Nodes name="all">
{nodes}
    </Nodes>
    <Elements type="penta6" name="e">
      <elem id="1">{",".join(str(c) for c in conn)}</elem>
    </Elements>
  </Mesh>
  <MeshDomains>
    <SolidDomain name="e" mat="m"/>
  </MeshDomains>
  <Boundary>
    <bc name="fix" type="zero displacement" node_set="n1">
      <x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof>
    </bc>
  </Boundary>
  <Control>
    <analysis>STATIC</analysis>
    <time_steps>1</time_steps>
    <step_size>1.0</step_size>
  </Control>
</febio_spec>
"""


def build_ns():
    return '<NodeSet name="n1">1,2,3</NodeSet>'


for name, conn in VARIANTS.items():
    d = TMP / name.split("_")[0]
    d.mkdir(exist_ok=True)
    f = d / "t.feb"
    txt = build(conn)
    txt = txt.replace('  <MeshDomains>',
                      f'    <NodeSet name="n1">1,2,3,4,5,6</NodeSet>\n'
                      f'  </MeshDomains>\n  <MeshDomains>')
    f.write_text(txt, encoding="utf-8")
    r = subprocess.run([str(FEBIO), "-i", f.name], cwd=d,
                       capture_output=True, text=True, timeout=120)
    out = (r.stdout or "") + (r.stderr or "")
    if "Negative jacobian" in out and "initialization" in out:
        verdict = "✗ 初始化负 Jacobian"
    elif "ERROR" in out:
        verdict = "? 其它错误: " + [l.strip() for l in out.split("\n")
                                   if "ERROR" in l][:1].__str__()
    else:
        verdict = "✓ 通过初始化"
    print(f"  {name:28s} {verdict}")
    if "✓" in verdict:
        # 看求解到哪一步
        for l in out.split("\n"):
            if "time steps completed" in l or "N O R M A L" in l:
                print(f"      {l.strip()[:80]}")
