"""看 FEBio 自己的 .log，找 Negative jacobian 附近是否有单元/域信息。"""
from pathlib import Path

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio_lig")
for nm in ("ankle_lig_WK.log",):
    p = W / nm
    if not p.exists():
        continue
    L = p.read_text(encoding="utf-8", errors="ignore").split("\n")
    print(f"=== {nm} ({len(L)} 行) ===")
    for i, l in enumerate(L):
        if "jacobian" in l.lower() or "negative" in l.lower():
            for k in range(max(0, i - 12), min(len(L), i + 6)):
                mark = ">>" if k == i else "  "
                print(f"{mark} {k+1:5d}| {L[k][:140]}")
            print()
            break
    # 也看尾部
    print("--- 尾部 25 行 ---")
    for k in range(max(0, len(L) - 25), len(L)):
        print(f"  {k+1:5d}| {L[k][:140]}")
