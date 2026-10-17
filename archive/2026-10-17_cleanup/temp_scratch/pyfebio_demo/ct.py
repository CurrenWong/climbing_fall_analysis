"""确认 pyvista 里线性四面体的 cell type 编号。"""

import pyvista as pv

for n in dir(pv.CellType):
    if "TETRA" in n:
        print(f"  {n:<28} = {int(getattr(pv.CellType, n))}")

print("\nVTK 约定: 10=VTK_QUADRATIC_TETRA(10节点), 12=VTK_TETRA(4节点)")
