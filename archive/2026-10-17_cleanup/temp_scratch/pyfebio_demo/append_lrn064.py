"""把 LRN-064 片段追加到权威 LEARNINGS.md（带撞号保护）。"""
import io
import os

SRC = r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\lrn064_snippet.md"
DST = r"D:\Project\climbing_fall_analysis\.learnings\LEARNINGS.md"
ID = "[LRN-20261005-064]"

s = io.open(SRC, encoding="utf-8").read()
before = os.path.getsize(DST)
t = io.open(DST, encoding="utf-8").read()
if ID in t:
    raise SystemExit(f"⚠ 撞号中止：{ID} 已存在，未写入")
io.open(DST, "a", encoding="utf-8").write(s)
after = os.path.getsize(DST)
t2 = io.open(DST, encoding="utf-8").read()
print(f"LEARNINGS.md {before} -> {after} bytes (+{after - before})")
print(f"{ID} 已写入? {ID in t2}   出现 {t2.count(ID)} 次")
