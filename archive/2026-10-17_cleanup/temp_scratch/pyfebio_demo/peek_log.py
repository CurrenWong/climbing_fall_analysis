"""检查 FEBio logfile 文本格式（步记录结构）。"""
import re
from pathlib import Path

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio_lig")
for f in ["disp_WK.txt", "reac_tibia_WK.txt", "reac_talus_WK.txt"]:
    p = W / f
    if not p.exists():
        print(f"(缺) {f}")
        continue
    L = p.read_text(encoding="utf-8", errors="replace").split("\n")
    recs = [(i, l) for i, l in enumerate(L) if l.startswith("*") or "Data" in l or "Step" in l]
    print(f"\n=== {f}  {len(L)} 行  记录头 {len(recs)} 个 ===")
    for i, l in recs[:6]:
        print(f"  {i:>5}|{l.strip()[:90]}")
    # 找最后一个记录的数据行
    last = [i for i, l in enumerate(L) if l.startswith("*Step")]
    if last:
        j = last[-1]
        print(f"  --- 最后一段 (行 {j+1}..) ---")
        for k in range(j, min(j + 8, len(L))):
            print(f"  {k:>5}|{L[k].strip()[:90]}")
        # 数据行数
        n = 0
        for k in range(j + 4, len(L)):
            if L[k].startswith("*") or not L[k].strip():
                break
            n += 1
        print(f"  --- 该段数据行数 {n} ---")
