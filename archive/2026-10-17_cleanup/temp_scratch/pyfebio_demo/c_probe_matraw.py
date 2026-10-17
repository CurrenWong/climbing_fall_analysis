"""查看问题 mid 的原始材料卡 —— 不猜字段位置，直接读原样。"""
import re
import zipfile
from pathlib import Path

ZIP = Path(r"D:\Project\climbing_fall_analysis\paper\AM50_V71_Occupant.zip")
MAIN = "AM50_V71_Occupant/THUMS_model/THUMS_AM50_V71_Occupant_202411.k"
WANT = {81000001, 81000601, 81000801, 81001000, 81001100, 81001200,
        81001300, 81000200, 81000300}

RE_KW = re.compile(r"^\*([A-Z_0-9]+)")
hits: dict[int, tuple[str, list[str]]] = {}
cur = None
blk: list[str] = []
blk_kw = None
buf = b""


def flush():
    global blk, blk_kw
    if blk and blk_kw and blk_kw.startswith("MAT_"):
        try:
            mid = int(blk[0][:10])
        except (ValueError, IndexError):
            blk = []
            return
        if mid in WANT and mid not in hits:
            hits[mid] = (blk_kw, list(blk))
    blk = []


with zipfile.ZipFile(ZIP) as z, z.open(MAIN) as f:
    while True:
        chunk = f.read(8 * 1024 * 1024)
        if not chunk:
            break
        buf += chunk
        txt = buf.decode("latin-1")
        lines = txt.split("\n")
        buf = lines.pop().encode("latin-1")
        for raw in lines:
            s = raw.rstrip("\r")
            st = s.strip()
            if not st or st[0] == "$":
                continue
            if st[0] == "*":
                flush()
                cur = st
                blk_kw = RE_KW.match(st.upper()).group(1) if RE_KW.match(st.upper()) else None
                continue
            if cur and cur.upper().startswith("*MAT"):
                blk.append(s)
    flush()

print(f"找到 {len(hits)}/{len(WANT)} 个 mid\n")
for mid in sorted(hits):
    kw, data = hits[mid]
    print(f"=== mid {mid}  ({kw}) ===")
    for i, l in enumerate(data[:4]):
        print(f"  card{i+1}| {l!r}")
        # 10 字符切片便于对照
        segs = [l[j*10:(j+1)*10].strip() for j in range(len(l)//10)]
        print(f"         fields: {segs}")
    print()

missing = WANT - set(hits)
if missing:
    print(f"⚠️ 未找到: {sorted(missing)}")
