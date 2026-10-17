"""pyfebio 官方 README 样例（单 hex8 压缩实验）。

用法:
    python demo.py            # 只生成 my_model.feb
    python demo.py --run      # 生成后再调 febio4 求解（需 FEBio 在 PATH）

来源: https://github.com/febiosoftware/pyfebio#example
"""

import sys
from pathlib import Path

import pyfebio

OUT = Path(__file__).parent / "my_model.feb"


def build() -> pyfebio.model.Model:
    my_model = pyfebio.model.Model()

    # 单个 hex8 单元的 8 个节点
    nodes_list = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0],
    ]
    elements_list = [[1, 2, 3, 4, 5, 6, 7, 8]]

    nodes = pyfebio.mesh.Nodes(name="nodes")
    for i, node in enumerate(nodes_list):
        nodes.add_node(pyfebio.mesh.Node(id=i + 1, text=",".join(map(str, node))))

    elements = pyfebio.mesh.Elements(name="box", type="hex8")
    for i, element in enumerate(elements_list):
        elements.add_element(
            pyfebio.mesh.Hex8Element(id=i + 1, text=",".join(map(str, element)))
        )

    my_model.mesh_.nodes.append(nodes)
    my_model.mesh_.elements.append(elements)

    bottom_nodes = [1, 2, 3, 4]
    top_nodes = [5, 6, 7, 8]
    my_model.mesh_.node_sets.append(
        pyfebio.mesh.NodeSet(name="top", text=",".join(map(str, top_nodes)))
    )
    my_model.mesh_.node_sets.append(
        pyfebio.mesh.NodeSet(name="bottom", text=",".join(map(str, bottom_nodes)))
    )

    # 材料: CoupledMooneyRivlin 软骨
    my_material = pyfebio.material.CoupledMooneyRivlin(
        id=1,
        name="cartilage",
        c1=pyfebio.material.MaterialParameter(text=10.0),
        c2=pyfebio.material.MaterialParameter(text=1.0),
        k=pyfebio.material.MaterialParameter(text=1000.0),
    )
    my_model.meshdomains_.add_solid_domain(
        pyfebio.meshdomains.SolidDomain(name="box", mat="cartilage")
    )
    my_model.material_.add_material(my_material)

    # 底面固定, 顶面 z 方向压缩 20%
    my_model.boundary_.add_bc(
        pyfebio.boundary.BCZeroDisplacement(
            node_set="bottom", x_dof=1, y_dof=1, z_dof=1
        )
    )
    displacement_value = pyfebio.boundary.Value(lc=1, text=-0.2)
    my_model.boundary_.add_bc(
        pyfebio.boundary.BCPrescribedDisplacement(
            node_set="top", dof="z", value=displacement_value
        )
    )

    # load curve 1: 时间 0->1, 比例 0->1
    my_model.loaddata_.add_load_curve(
        pyfebio.loaddata.LoadCurve(
            id=1, points=pyfebio.loaddata.CurvePoints(points=["0.0,0.0", "1.0,1.0"])
        )
    )
    return my_model


def main() -> int:
    model = build()
    model.save(OUT)
    print(f"[ok] 已生成 {OUT}  ({OUT.stat().st_size} bytes)")

    if "--run" not in sys.argv:
        print("[skip] 未指定 --run, 跳过 FEBio 求解")
        return 0

    print("[run] 调用 febio4 ...")
    rc = pyfebio.model.run_model(OUT, silent=True)
    print(f"[run] febio4 returncode = {rc}  (非 0 通常表示 FEBio 未安装或求解失败)")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
