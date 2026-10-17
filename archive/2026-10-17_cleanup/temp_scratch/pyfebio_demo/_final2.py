import io
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis")
for f in ("scripts/ankle_fe/thums_materials.py", "scripts/ankle_fe/cartilage.py",
          "scripts/ankle_fe/build_feb.py", "temp/thums/bone_materials.json",
          "temp/thums/cartilage.json", "temp/thums/cartilage_mesh.npz",
          "temp/thums/febio/ankle_tibiotalar.feb", "temp/thums/febio/febio.log",
          "FE_PIPELINE_TODO.md", "THUMS_EXTRACT_TODO.md",
          "docs/踝关节FE建模交接.md", ".learnings/LEARNINGS.md"):
    p = P / f
    print(f"  {'OK ' if p.exists() else 'MISS'} {p.stat().st_size if p.exists() else 0:>10,}  {f}")

t = io.open(P / "docs/踝关节FE建模交接.md", encoding="utf-8").read()
need = ["## 13.", "13.1", "13.2", "13.3", "13.4", "12.5", "12.6"]
print("\n交接文档缺失:", [n for n in need if n not in t] or "无 OK",
      f"({len(t.encode('utf-8'))} bytes)")

lrn = io.open(P / ".learnings/LEARNINGS.md", encoding="utf-8").read()
print("LRN-014:", "LRN-20261003-014" in lrn, " LRN-015:", "LRN-20261003-015" in lrn,
      f"({len(lrn.encode('utf-8'))} bytes)")
