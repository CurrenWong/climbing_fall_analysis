"""让 pyfebio 自己生成含 TractionLoad / Contact 的 .feb，读回 XML 学确切语法。"""
from pathlib import Path
import pyfebio

print("Model 字段:", list(pyfebio.model.Model.model_fields))
m = pyfebio.model.Model()
print("\n各 section 属性:", [n for n in dir(m) if not n.startswith('_')])

# ---- 看 Literal 允许的 type 值 ----
from typing import get_args
for cls in (pyfebio.loads.TractionLoad, pyfebio.loads.NodalForce,
            pyfebio.loads.PressureLoad, pyfebio.loads.NodalLoad):
    f = cls.model_fields.get("type")
    print(f"{cls.__name__:16s} tag={getattr(cls,'xml_tag',None)} "
          f"type Literal={get_args(f.annotation) if f else None}")

print("\nLoads 字段细节:")
for n, f in pyfebio.loads.Loads.model_fields.items():
    print(f"  {n}: annotation={f.annotation}")

print("\nContact 类:")
import pyfebio.contact as C
for n in dir(C):
    o = getattr(C, n)
    if isinstance(o, type):
        t = getattr(o, "xml_tag", None)
        if t:
            print(f"  {n:26s} tag={t} fields={list(o.model_fields)[:9]}")
