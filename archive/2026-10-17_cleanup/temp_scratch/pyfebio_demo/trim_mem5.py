"""最后 107 B：紧凑化 gmsh 与 P1 两行。"""
import io

MEM = r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md"
m = io.open(MEM, encoding="utf-8").read()
n0 = len(m.encode("utf-8"))
REPS = [
    ("- **`uv pip install gmsh` ≠ 官方 SDK**（occ 绑定残缺 ⇒ STL→实体不通）；"
     "官方在 `D:\\Program\\gmsh-sdk\\`，DLL 在 `lib/`。**同版本号可能是两份二进制；"
     "判据=实测 API 可调用性**",
     "- **`uv pip install gmsh` ≠ 官方 SDK**（occ 残缺 ⇒ STL→实体不通）；官方在 "
     "`D:\\Program\\gmsh-sdk\\`（DLL 在 `lib/`）。**同版本号可能是两份二进制，判据=实测 API**"),
    ("- **P1 硬结果**：抱石场景**一例骨折都报不出来**（1000 MC，最高利用率 0.725）"
     "⇒ 损伤分布**空**（无旋转 32.6%→可产类型 19.3%）",
     "- **P1 硬结果**：抱石场景**一例骨折都报不出**（1000 MC，最高利用率 0.725）"
     "⇒ 损伤分布**空**（无旋转 32.6%→19.3%）"),
]
for old, new in REPS:
    if old in m:
        m = m.replace(old, new)
    else:
        print(f"⚠ 未命中：{old[:40]}")
io.open(MEM, "w", encoding="utf-8").write(m)
n1 = len(m.encode("utf-8"))
print(f"MEMORY.md {n0} -> {n1} B  {'✅ ≤10240' if n1 <= 10240 else f'❌ 超 {n1-10240}'}")
