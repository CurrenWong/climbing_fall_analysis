"""全量扫描 THUMS 主 .k（230MB），统计关键字 + 抽取下肢骨部件名 + 判断有无节点/实体单元。
"""
import re
import zipfile
from collections import Counter

ZIP = r"D:\Project\climbing_fall_analysis\paper\AM50_V71_Occupant.zip"
INNER = "AM50_V71_Occupant/THUMS_model/THUMS_AM50_V71_Occupant_202411.k"

kw = Counter()
parts = []            # (title, id)
pending = None        # *PART 标题下一行是 part 名
in_part = False
nbytes = 0

BONE_HINT = re.compile(
    r"tibia|fibula|femur|patella|talus|calcaneus|calcn|navicular|"
    r"cuboid|metatars|cuneiform|phal|hallux|foot",
    re.I,
)

with zipfile.ZipFile(ZIP) as z, z.open(INNER) as f:
    buf = b""
    while True:
        chunk = f.read(8 * 1024 * 1024)
        if not chunk:
            break
        nbytes += len(chunk)
        buf += chunk
        txt = buf.decode("latin-1")
        lines = txt.split("\n")
        buf = lines.pop().encode("latin-1")

        for ln in lines:
            s = ln.strip()
            if not s:
                continue
            if s.startswith("$"):
                name = s[1:].strip().split()[0] if len(s) > 1 else ""
                kw["$" + name.upper()] += 1
                if s.startswith("$HMNAME COMPS"):
                    rest = s[len("$HMNAME COMPS"):].strip()
                    m = re.match(r"^(\d+)\s*(.*)$", rest)
                    if m:
                        parts.append((m.group(2).strip(), m.group(1)))
                continue
            if s.startswith("*"):
                name = s.split()[0].upper()
                kw[name] += 1
        # 保留最后半行
        del lines

print(f"全文件扫描 {nbytes/1024/1024:.1f} MB")
print("\n=== 全部关键字（按出现次数）===")
for k, v in kw.most_common():
    if k.startswith("*"):
        print(f"  {k:38s} {v:7d}")

print("\n=== 下肢骨相关部件（去重，按编号）===")
bone_parts = {}
for title, pid in parts:
    if BONE_HINT.search(title):
        bone_parts[pid] = title
for pid in sorted(bone_parts):
    print(f"  {pid}  {bone_parts[pid]}")

print(f"\n总部件 {len(parts)}，其中骨相关 {len(bone_parts)}")
print("\n=== 有无网格数据 ===")
print(f"  *NODE          : {kw.get('*NODE', 0)}")
print(f"  *ELEMENT_SOLID : {kw.get('*ELEMENT_SOLID', 0)}")
print(f"  *ELEMENT_SHELL : {kw.get('*ELEMENT_SHELL', 0)}")
print(f"  *INCLUDE       : {kw.get('*INCLUDE', 0)}")
