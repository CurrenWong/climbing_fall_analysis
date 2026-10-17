"""定位: pyfebio 0.3.0 反序列化 Literal[0,1] 字段的 bug。"""

from pydantic_xml import BaseXmlModel, attr, element
from typing import Literal


class A(BaseXmlModel, tag="bc"):
    type: Literal["zero displacement"] = attr(default="zero displacement", frozen=True)
    node_set: str = attr()
    x_dof: Literal[0, 1] = element(default=0)


XML = '<bc type="zero displacement" node_set="bottom"><x_dof>1</x_dof></bc>'

# 1) 序列化: pydantic 写 int -> "1"
a = A(node_set="bottom", x_dof=1)
print("serialize ->", a.to_xml())

# 2) 反序列化: lxml 给的是 str "1"
try:
    A.from_xml(XML)
    print("from_xml OK")
except Exception as e:
    print("from_xml FAILED:", str(e).splitlines()[1].strip())

# 3) 证明根因: Literal 收不了 str
from pydantic import BaseModel, ValidationError


class M(BaseModel):
    v: Literal[0, 1]


for raw in (1, "1", 0, "0"):
    try:
        M(v=raw)
        print(f"  Literal[0,1] <- {raw!r}: OK")
    except ValidationError:
        print(f"  Literal[0,1] <- {raw!r}: FAIL")

# 4) 对照: int 字段能收 str "1"
class N(BaseModel):
    v: int


print("  int <- '1':", N(v="1").v)
