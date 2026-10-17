"""生成含 TractionLoad + Contact 的 .feb，读 XML 学确切语法（v2）。"""
from pathlib import Path
import pyfebio
import pyfebio.contact as C

OUT = Path(__file__).parent / "probe_pyfebio_load.feb"

print("=== contact 模块成员 ===")
for n in dir(C):
    o = getattr(C, n)
    if isinstance(o, type) and n[0].isupper():
        print(f"  {n:30s} tag={getattr(o,'xml_tag',None)} fields={list(getattr(o,'model_fields',{}))[:10]}")

m = pyfebio.model.Model()
nodes = pyfebio.mesh.Nodes(name="nodes")
for i, p in enumerate([[0,0,0],[1,0,0],[1,1,0],[0,1,0],
                       [0,0,1],[1,0,1],[1,1,1],[0,1,1]], 1):
    nodes.add_node(pyfebio.mesh.Node(id=i, text=",".join(map(str, p))))
el = pyfebio.mesh.Elements(name="box", type="hex8")
el.add_element(pyfebio.mesh.Hex8Element(id=1, text="1,2,3,4,5,6,7,8"))
surf = pyfebio.mesh.Elements(name="topsurf", type="quad4")
surf.add_element(pyfebio.mesh.Quad4Element(id=1, text="5,6,7,8"))
m.mesh_.nodes.append(nodes)
m.mesh_.elements.extend([el, surf])
m.mesh_.node_sets.append(pyfebio.mesh.NodeSet(name="bottom", text="1,2,3,4"))
m.material_.add_material(pyfebio.material.IsotropicElastic(
    id=1, name="m",
    density=pyfebio.material.MaterialParameter(text=1e-6),
    E=pyfebio.material.MaterialParameter(text=1000),
    v=pyfebio.material.MaterialParameter(text=0.3)))
m.meshdomains_.add_solid_domain(pyfebio.meshdomains.SolidDomain(name="box", mat="m"))
m.boundary_.add_bc(pyfebio.boundary.BCZeroDisplacement(
    node_set="bottom", x_dof=1, y_dof=1, z_dof=1))

print("\n=== TractionLoad 构造 ===")
for kw in ({"traction": "0,0,-1"},):
    try:
        tl = pyfebio.loads.TractionLoad(name="pull", type="traction",
                                       surface="topsurf", **kw)
        m.loads_.all_surface_loads.append(tl)
        print(f"  [ok] {kw}")
    except Exception as exc:
        print(f"  [fail] {kw}: {type(exc).__name__}: {exc}")

print("\n=== Contact 构造 ===")
try:
    m.contact_.all_contact_interfaces.append(
        C.SlidingElastic(name="tib", surface_pair="topsurf,topsurf",
                         laugon="PENALTY", penalty=1.0, tolerance=0.1))
    print("  [ok] SlidingElastic")
except Exception as exc:
    print("  [fail]:", type(exc).__name__, exc)

try:
    m.save(OUT)
    xml = OUT.read_text(encoding="latin-1")
    print("\n=== 生成的 XML 尾部（Loads / Boundary / Contact）===")
    i = xml.find("<Loads")
    print(xml[i:])
except Exception as exc:
    print("save 失败:", type(exc).__name__, exc)
