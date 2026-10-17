"""看 P2① subagent 留下什么（按修改时间看最近产物）。"""
import os
import time
from pathlib import Path

R = Path(r"D:\Project\climbing_fall_analysis")
now = time.time()
CUT = 3600          # 最近 1 小时

DIRS = [R / "temp/pyfebio_demo", R / "temp/thums", R / "temp/thums/cart_m1_fix",
        R / "scripts/ankle_fe", R / "temp/thums/febio_lig"]
print("=" * 96)
print("最近 60 分钟内被修改的文件")
print("=" * 96)
rows = []
for d in DIRS:
    if not d.is_dir():
        print(f"  ✗ 不存在: {d}")
        continue
    for f in d.iterdir():
        if f.is_file():
            m = f.stat().st_mtime
            if now - m < CUT:
                rows.append((m, f))
rows.sort()
for m, f in rows:
    age = now - m
    print(f"  {time.strftime('%H:%M:%S', time.localtime(m))}  "
          f"({age/60:5.1f} 分钟前)  {f.stat().st_size:>10,} B  "
          f"{f.relative_to(R)}")

print("\n" + "=" * 96)
print("cart_m1_fix 目录")
print("=" * 96)
d = R / "temp/thums/cart_m1_fix"
if d.is_dir():
    for f in sorted(d.iterdir()):
        print(f"  {f.name:<44} {f.stat().st_size:>10,} B")
else:
    print('  ✗ 未创建 ⇒ subagent 还没走到「写新网格」这步')

print("\n" + "=" * 96)
print("cart_patch.py 是否被改过")
print("=" * 96)
p = R / "scripts/ankle_fe/cart_patch.py"
m = p.stat().st_mtime
print(f"  {p.name}: {p.stat().st_size:,} B, 修改于 "
      f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(m))} "
      f"({(now-m)/60:.1f} 分钟前)")
print(f"  ⇒ {'**被改过**' if now-m < CUT else '未被改过（仍是原版）'}")
