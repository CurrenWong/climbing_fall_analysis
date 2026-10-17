"""交接完整性最终验证（用正确的项目根）。"""
import re
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")

t = (ROOT / "FE_PIPELINE_TODO.md").read_text(encoding="utf-8")
head = t[t.index("# 🧭 接手入口"):t.index("## C. 解出 THUMS")]

paths = sorted(set(re.findall(r"`([\w\-./\\]+\.(?:py|md|json|npz|zip|k|fe?b))`",
                              head)))
miss = []
for p in paths:
    if (ROOT / p.replace("\\", "/")).exists():
        continue
    if (WS / p).exists():
        continue
    miss.append(p)
print(f"=== 入口节路径：{len(paths)} 个，{len(miss)} 个解析不了 ===")
for p in paths:
    print(f"  {'OK  ' if p not in miss else 'MISS'} {p}")

print("\n=== 开关表不再有裸文件格 ===")
print("  ", "| `build_feb.py` |" not in t and "| `cartilage.py` |" not in t)

print("\n=== 开关表渲染完整性（每行 3 格）===")
tbl = t[t.index("## 环境开关"):t.index("## 关键文件")]
for l in tbl.split("\n"):
    if l.startswith("| `"):
        n = l.count(" | ") + 1
        flag = "OK" if n == 3 else f"❌ {n} 格"
        print(f"  {flag}  {l[:70]}")

print("\n=== 全量收尾检查 ===")
checks = {
    "FE_PIPELINE 有接手入口": "# 🧭 接手入口" in t,
    "入口标注已实跑验证": "本节命令已实跑验证" in t,
    "过时状态已清除": "A ⚠️ 进行中" not in t and "卡在 Jacobian" not in t,
    "ANKLE_FE 标为已取代": "IITD 路线已停止" in
        (ROOT / "ANKLE_FE_TODO.md").read_text(encoding="utf-8"),
    "交接文档互指入口": "先读 `FE_PIPELINE_TODO.md` 的「🧭 接手入口」节" in
        (ROOT / "docs/踝关节FE建模交接.md").read_text(encoding="utf-8"),
    "MEMORY ≤ 10KB": len((WS / "MEMORY.md").read_text(
        encoding="utf-8").encode("utf-8")) <= 10240,
    "LRN-016~020 齐": all(f"LRN-20261003-0{n}" in
        (ROOT / ".learnings/LEARNINGS.md").read_text(encoding="utf-8")
        for n in (16, 17, 18, 19, 20)),
}
for k, v in checks.items():
    print(f"  {'✅' if v else '❌'} {k}")
print(f"\nFE_PIPELINE_TODO.md {len(t.encode('utf-8')):,} bytes")
