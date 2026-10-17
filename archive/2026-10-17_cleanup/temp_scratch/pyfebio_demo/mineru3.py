"""完整（不截断）导出关键段落。"""
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis\paper"
         r"\MinerU_markdown_人体足踝部有限元模型的建立及有效性分析_吴恺"
         r"_2106639040860729344.md")
L = P.read_text(encoding="utf-8", errors="ignore").splitlines()

KEY = [94, 114, 135, 137, 139, 15, 17]
out = []
for n in KEY:
    if 1 <= n <= len(L):
        out.append(f"\n{'='*76}\n行 {n}（完整）\n{'='*76}")
        out.append(L[n - 1])
out.append(f"\n\n{'='*76}\n全文不截断（共 {len(L)} 行）\n{'='*76}")
for i, l in enumerate(L, 1):
    out.append(f"{i:4d}| {l}")
Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\_mu_raw.txt"
     ).write_text("\n".join(out), encoding="utf-8")
print("wrote", len("\n".join(out)))
