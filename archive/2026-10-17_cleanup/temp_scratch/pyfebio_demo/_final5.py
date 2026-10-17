"""第二轮收尾：语法/产物/开关 全面自检。"""
import ast
import os
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")

print("=== 脚本语法 ===")
for f in ("scripts/ankle_fe/cartilage.py", "scripts/ankle_fe/build_feb.py",
          "scripts/ankle_fe/thums_materials.py", "scripts/ankle_fe/thums_extract.py"):
    try:
        ast.parse((P / f).read_text(encoding="utf-8"))
        print(f"  OK   {f}")
    except SyntaxError as e:
        print(f"  FAIL {f}: {e}")

print("\n=== 产物 ===")
for f in ("temp/thums/cartilage.json", "temp/thums/cartilage_mesh.npz",
          "temp/thums/febio/ankle_tibiotalar.feb",
          "temp/thums/febio/ankle_tibiotalar.log",
          "temp/thums/febio/ankle_tibiotalar.xplt",
          "FE_PIPELINE_TODO.md", "docs/踝关节FE建模交接.md",
          ".learnings/LEARNINGS.md",
          "temp/pyfebio_demo/diag_facets.py"):
    p = P / f
    print(f"  {'OK ' if p.exists() else 'MISS'} {p.stat().st_size if p.exists() else 0:>10,}  {f}")

print("\n=== 环境开关默认值 ===")
t = (P / "scripts/ankle_fe/cartilage.py").read_text(encoding="utf-8")
b = (P / "scripts/ankle_fe/build_feb.py").read_text(encoding="utf-8")
for k, src, name in (("CART_QUALITY", t, "cartilage"), ("CAP_SAFETY", t, "cartilage"),
                     ("FULL_ANKLE", b, "build_feb"), ("LAUGON", b, "build_feb"),
                     ("NODE_RELOC", b, "build_feb"), ("TWO_PASS", b, "build_feb"),
                     ("PENALTY", b, "build_feb"), ("TIME_STEPS", b, "build_feb"),
                     ("DRIVER", b, "build_feb"), ("DISP_MM", b, "build_feb")):
    line = next((l.strip() for l in src.split("\n") if l.strip().startswith(k)), None)
    print(f"  {name:9s} {line[:78] if line else 'MISSING'}")

print("\n=== 文档 ===")
m = (WS / "memory/2026-10-03/ankle-fe-modeling.md").read_text(encoding="utf-8")
l = (WS / ".learnings/LEARNINGS.md").read_text(encoding="utf-8")
mm = (WS / "MEMORY.md").read_text(encoding="utf-8")
print(f"  记忆 §九/§十: {'## 九、' in m} / {'## 十、' in m}"
      f"  ({len(m.encode('utf-8')):,} bytes)")
print(f"  LRN-016..020: "
      f"{all(f'LRN-20261003-0{n}' in l for n in (16, 17, 18, 19, 20))}")
print(f"  MEMORY.md {len(mm.encode('utf-8')):,} bytes（自设上限 10240，⚠️ 待整理）")
