"""追加 LRN-070（先查重号，避免静默跳过）。"""
import io
import re

LEARN = r"D:\Project\climbing_fall_analysis\.learnings\LEARNINGS.md"
SNIP = r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\lrn070_snippet.md"
txt = io.open(LEARN, encoding="utf-8").read()
new_id = "LRN-20261005-070"
hits = re.findall(re.escape(new_id), txt)
print(f"查重：{new_id} 已出现 {len(hits)} 次")
if hits:
    print("★ 撞号 ⇒ 不写，需换号")
else:
    snip = io.open(SNIP, encoding="utf-8").read()
    with io.open(LEARN, "a", encoding="utf-8") as fh:
        fh.write(snip)
    print(f"已追加 {len(snip)} 字符")
    after = io.open(LEARN, encoding="utf-8").read()
    print(f"写入后：{new_id} 出现 {after.count(new_id)} 次；文件 {len(after)} 字符")
