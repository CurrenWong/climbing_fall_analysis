"""把插入块移到「🧭 接手入口」标题之后（更自然）。"""
from pathlib import Path

FP = Path(r"D:\Project\climbing_fall_analysis\FE_PIPELINE_TODO.md")
L = FP.read_text(encoding="utf-8").split("\n")

bi = next(i for i, l in enumerate(L) if l.startswith("> 🆕 **2026-10-04 论文重解析"))
blk = [L[bi], L[bi + 1]]          # 提示行 + 空行
del L[bi:bi + 2]

hi = next(i for i, l in enumerate(L) if "接手入口" in l or "\U0001f9ed" in l)
L[hi + 1:hi + 1] = [""] + blk
FP.write_text("\n".join(L), encoding="utf-8")

print(f"{FP.stat().st_size:,} bytes")
for i in range(max(0, hi - 3), hi + 6):
    print(f"{i+1:4d}| {L[i][:118]}")
