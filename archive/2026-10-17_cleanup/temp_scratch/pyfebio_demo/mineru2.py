"""dump MinerU 解析的关键段落：图4 / 结果 / 表格。"""
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis\paper"
         r"\MinerU_markdown_人体足踝部有限元模型的建立及有效性分析_吴恺"
         r"_2106639040860729344.md")
t = P.read_text(encoding="utf-8", errors="ignore")
L = t.splitlines()

out = []
for i, l in enumerate(L, 1):
    if "图4" in l or "表1" in l:
        i0, i1 = max(0, i - 4), min(len(L), i + 5)
        out.append(f"\n{'='*76}\n[图4/表1] 行 {i}\n{'='*76}")
        for k in range(i0, i1):
            out.append(f"{k+1:4d}| {L[k][:190]}")

out.append(f"\n\n{'='*76}\n全文（310 行）\n{'='*76}")
for i, l in enumerate(L, 1):
    out.append(f"{i:4d}| {l[:190]}")

Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\_mu_full.txt"
     ).write_text("\n".join(out), encoding="utf-8")
print("wrote", len("\n".join(out)))
