"""探测 THUMS .k 的字段格式：单次流式扫描，抓三类块的原始样本 + 块位置。

目的：LS-DYNA 关键字卡有「逗号分隔」与「固定宽度(I8/8F)」两种写法，
猜错会使 eid/pid/nid 全部读错 —— 必须先看原样。
"""
import zipfile

ZIP = r"D:\Project\climbing_fall_analysis\paper\AM50_V71_Occupant.zip"
INNER = "AM50_V71_Occupant/THUMS_model/THUMS_AM50_V71_Occupant_202411.k"

WANT = {
    "*NODE": 14,
    "*ELEMENT_SOLID": 14,
    "*ELEMENT_SHELL": 14,
    "*PART": 14,
    "*SECTION_SOLID": 6,
}

samples = {k: [] for k in WANT}
offsets = {}          # kw -> (byte offset of block start, count)
counts = {}
block_start = {}
cur = None
ln_no = 0
pos = 0

with zipfile.ZipFile(ZIP) as z, z.open(INNER) as f:
    buf = b""
    while True:
        chunk = f.read(8 * 1024 * 1024)
        if not chunk:
            break
        buf += chunk
        txt = buf.decode("latin-1")
        lines = txt.split("\n")
        tail = lines.pop()
        start_pos = pos
        pos += len(chunk)
        buf = tail.encode("latin-1")

        for i, ln in enumerate(lines):
            ln_no += 1
            s = ln.rstrip("\r")
            st = s.strip()
            if st.startswith("*"):
                kw = st.split()[0].upper()
                counts[kw] = counts.get(kw, 0) + 1
                if kw in WANT:
                    cur = kw
                    block_start.setdefault(kw, (ln_no, start_pos))
                    continue
                else:
                    cur = None
                    continue
            if cur and st and not st.startswith("$"):
                if len(samples[cur]) < WANT[cur]:
                    samples[cur].append((ln_no, repr(s)))

print("=== 块首次出现位置 ===")
for k in WANT:
    print(f"  {k:18s} first at line {block_start.get(k, ('-','-'))[0]}")

print("\n=== 各关键字出现次数（已知）===")
for k in ("*NODE", "*ELEMENT_SOLID", "*ELEMENT_SHELL", "*PART", "*SECTION_SOLID"):
    print(f"  {k:18s} {counts.get(k, 0)}")

for k in WANT:
    print(f"\n=== {k} 前 {len(samples[k])} 行数据（原样）===")
    for ln_no, raw in samples[k]:
        print(f"  L{ln_no:<8d} {raw}")
