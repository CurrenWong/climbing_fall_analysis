"""往返验证: from_xml 读回刚生成的 .feb, 确认能被自己解析。"""

from pathlib import Path

import pyfebio

feb = Path(__file__).parent / "my_model.feb"
# 注意: .feb 带 <?xml version='1.0' encoding='ISO-8859-1'?> 声明,
# lxml 拒收 "带编码声明的 str", 必须传 bytes;
# 而 pydantic 的 parse_file 是 JSON 解析, 对 XML 不适用。
m = pyfebio.model.Model.from_xml(feb.read_bytes())

print("parse ok")
print("  material   :", [(x.name, x.type) for x in m.material_.materials])
print("  nodes      :", len(m.mesh_.nodes[0].nodes))
print("  elements   :", m.mesh_.elements[0].type, len(m.mesh_.elements[0].elements))
print("  node_sets  :", [ns.name for ns in m.mesh_.node_sets])
print("  bcs        :", [b.type for b in m.boundary_.bcs])
print("  load curves:", [c.id for c in m.loaddata_.load_curves])
