"""在 FEBio 自己的 .log 里找 Negative jacobian 的上下文（可能有单元号）。"""
from pathlib import Path

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio_lig")
p = W / "ankle_lig_WK.log"
L = p.read_text(encoding="utf-8", errors="ignore").split("\n")
print(f"{p.name}: {len(L)} 行")

hits = [i for i, l in enumerate(L)
        if "jacobian" in l.lower() or "negative" in l.lower()]
print(f"命中 {len(hits)} 行")
for i in hits[:20]:
    print(f"\n--- 行 {i+1} 上下文 ---")
    for k in range(max(0, i - 6), min(len(L), i + 8)):
        print(f"  {k+1:5d}| {L[k][:150]}")
