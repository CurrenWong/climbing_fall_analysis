"""读 cart_m1_fix 的验收报告 + 与基线对照。"""
import io
import re

R = r"D:\Project\climbing_fall_analysis"
for tag, f in (("基线 cart_m1", "_reg.txt"), ("新 cart_m1_fix", "_accfix.txt")):
    p = f"{R}/temp/pyfebio_demo/{f}"
    t = io.open(p, encoding="utf-8", errors="replace").read()
    print("=" * 92)
    print(f"=== {tag}   ({len(t)} 字符)")
    print("=" * 92)
    for l in t.splitlines():
        if re.search(r"求解|末态|接触面|接触\s*:|Loads|合计|单元：|写出|负 ?J|NORMAL|ERROR"
                     r"|合力|u\|max|刚体", l):
            print("  " + l.strip()[:136])
    print()
