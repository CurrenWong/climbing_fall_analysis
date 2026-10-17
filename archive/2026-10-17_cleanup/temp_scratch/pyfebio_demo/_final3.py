import io
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis")
MEM = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md")

print("=== 文件 ===")
for f in ("scripts/ankle_fe/thums_materials.py", "scripts/ankle_fe/cartilage.py",
          "scripts/ankle_fe/build_feb.py", "temp/thums/bone_materials.json",
          "temp/thums/cartilage.json", "temp/thums/cartilage_mesh.npz",
          "temp/thums/febio/ankle_tibiotalar.feb", "temp/thums/febio/ankle_tibiotalar.xplt",
          "FE_PIPELINE_TODO.md", "docs/踝关节FE建模交接.md", ".learnings/LEARNINGS.md"):
    p = P / f
    print(f"  {'OK ' if p.exists() else 'MISS'} {p.stat().st_size if p.exists() else 0:>10,}  {f}")

t = io.open(P / "FE_PIPELINE_TODO.md", encoding="utf-8").read()
print(f"\nTODO {len(t.encode('utf-8'))} bytes；章节:",
      [l.strip() for l in t.split("\n") if l.startswith("## ")])
d = io.open(P / "docs/踝关节FE建模交接.md", encoding="utf-8").read()
print(f"交接文档 {len(d.encode('utf-8'))} bytes")
for k in ("13.1", "13.2", "13.3", "13.4", "13.5", "13.6", "13.7"):
    if k not in d:
        print(f"  缺失 {k}")
l = io.open(P / ".learnings/LEARNINGS.md", encoding="utf-8").read()
print("LRN 016/017/018:",
      all(f"LRN-20261003-0{n}" in l for n in (16, 17, 18)))

m = MEM.read_text(encoding="utf-8")
print(f"\nMEMORY.md {len(m.encode('utf-8'))} bytes (上限 10240): "
      f"{'OK' if len(m.encode('utf-8')) <= 10240 else '⚠️ 超标'}")
print("gotcha 行数:", sum(1 for x in m.split("\n") if x.startswith("| 2") or
                          x.startswith("| 1")))
