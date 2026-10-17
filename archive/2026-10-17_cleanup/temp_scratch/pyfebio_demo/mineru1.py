"""看 MinerU 重解析的论文 markdown：找图4/韧带位移的原始数据。"""
import re
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis\paper"
         r"\MinerU_markdown_人体足踝部有限元模型的建立及有效性分析_吴恺"
         r"_2106639040860729344.md")
print(f"存在: {P.exists()}")
if not P.exists():
    d = P.parent
    for q in sorted(d.glob("*.md")):
        print(f"  {q.stat().st_size:>10,}  {q.name}")
    raise SystemExit(0)

t = P.read_text(encoding="utf-8", errors="ignore")
print(f"字符 {len(t):,}")
print(f"行数 {len(t.splitlines()):,}")

print("\n=== 关键数字是否出现 ===")
for kw in ("5.24", "5.04", "2.05", "1.92", "240", "260", "7300", "600",
           "3.74", "图 4", "图4", "图 3", "表 1", "表1"):
    n = len(re.findall(re.escape(kw), t))
    print(f"  {kw!r:10s} {n} 处")

print("\n=== 含「韧带」的行（前 40）===")
c = 0
for i, l in enumerate(t.splitlines(), 1):
    if "韧带" in l:
        print(f"  {i:4d}| {l.strip()[:150]}")
        c += 1
        if c >= 40:
            break

print("\n=== 图片/表格标记 ===")
for m in re.finditer(r"!\[[^\]]*\]\([^)]*\)|\|.*\|", t):
    s = m.group(0)
    if len(s) < 200:
        print("   ", s[:160])
