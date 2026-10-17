"""gmsh .geo 脚本法: Merge 离散 STL -> Surface Loop -> Volume -> 四面体网格。

依据 gmsh 官方 FAQ:
  "You can import a surface mesh in any one of the supported mesh file formats,
   define a volume, and mesh it.  See examples/simple_geo/sphere-discrete.geo"

关键差异: .geo 脚本里的 Merge() 会为离散实体建参数化几何,
而 Python API 的 gmsh.merge() 只建 "Discrete surface"(无 CAD 描述),
导致 geo.addSurfaceLoop 拿不到可闭合的环 -> 无法 addVolume。
用 .geo 脚本让 gmsh 自己走完整流程。
"""

import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
GEOM = ROOT / "temp" / "ankle_geom"
OUT = ROOT / "temp" / "ankle_mesh"
OUT.mkdir(parents=True, exist_ok=True)

GEO_TMPL = """// {name}
Merge "{stl}";
// Merge 对离散面建了参数化几何, 之后可枚举曲线/曲面
Curve Loop(1) = {{1}};
Surface {{1}} = {{1}};
Volume(1) = {{1}};
"""

# 需要先探测 Merge 之后实体 tag 的分配, 所以分两步
PROBE_TMPL = """
Merge "{stl}";
Print "ENTITIES";
Curve Loop(9999) = {{1}};
"""


def probe(stl: Path) -> dict:
    """跑一次最小 .geo, 打印 gmsh 认为的实体结构。"""
    geo = OUT / "_probe.geo"
    geo.write_text(f'Merge "{stl.as_posix()}";\n', encoding="utf-8")
    p = subprocess.run(
        ["gmsh", str(geo), "-parse_and_exit", "-v", "3"],
        capture_output=True,
        text=True,
        timeout=300,
    )
    return {"stdout": p.stdout, "stderr": p.stderr, "code": p.returncode}


if __name__ == "__main__":
    stl = GEOM / (sys.argv[1] if len(sys.argv) > 1 else "talus_r.stl")
    print(f"STL: {stl}")

    # 先看 gmsh 命令行版在不在
    exe = subprocess.run(["where", "gmsh"], capture_output=True, text=True)
    print(f"gmsh.exe: {exe.stdout.strip() or '<未安装>'}")
    if not exe.stdout.strip():
        sys.exit(1)

    r = probe(stl)
    print(f"\nrc={r['code']}")
    print("--- stdout ---")
    print(r["stdout"][-3000:])
    if r["stderr"]:
        print("--- stderr ---")
        print(r["stderr"][-2000:])
