"""看 target_cmp.py / final_600.py / gap_area.py 是不是"双靶对比"工具。"""
import io
import os

R = r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo"
for f in ("target_cmp.py", "final_600.py", "gap_area.py"):
    p = os.path.join(R, f)
    if not os.path.exists(p):
        print(f"✗ {f}")
        continue
    t = io.open(p, encoding="utf-8", errors="ignore").read()
    ls = t.splitlines()
    print("=" * 92)
    print(f"=== {f}   {len(t.encode('utf-8'))} B / {len(ls)} 行")
    print("=" * 92)
    for l in ls[:26]:
        print("  " + l[:136])
    print("  ...")
    for l in ls[-12:]:
        print("  " + l[:136])
