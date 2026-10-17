"""查 THUMS 是否已有独立软骨部件（之前搜索模式漏了 cart/chondr 等关键词）。"""
import re
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

ZIP = Path(r"D:\Project\climbing_fall_analysis\paper\AM50_V71_Occupant.zip")
MAIN = "AM50_V71_Occupant/THUMS_model/THUMS_AM50_V71_Occupant_202411.k"
PAT = re.compile(r"cart|chondr|menisc|synov|labr|disc|soft|tissue", re.I)
RE_KW = re.compile(r"^\*([A-Z_0-9]+)")

parts: dict[int, str] = {}
pending = None
elem = defaultdict(Counter)     # kw -> Counter(pid)
cur = None
buf = ""
FIELD = 8


def sl(s, w):
    v = []
    for i in range(0, len(s) - len(s) % w, w):
        f = s[i:i + w].strip()
        if not f:
            continue
        try:
            v.append(int(f))
        except ValueError:
            return None
    return v


with zipfile.ZipFile(ZIP) as z, z.open(MAIN) as f:
    while True:
        chunk = f.read(8 * 1024 * 1024)
        if not chunk:
            break
        buf += chunk.decode("latin-1")
        lines = buf.split("\n")
        buf = lines.pop()
        for raw in lines:
            s = raw.rstrip("\r")
            st = s.strip()
            if not st or st[0] == "$":
                continue
            if st[0] == "*":
                cur = RE_KW.match(st.upper()).group(1)
                if cur == "PART":
                    pending = None
                continue
            if cur == "PART":
                if pending is None:
                    pending = st
                else:
                    v = sl(st, 10)
                    if v:
                        parts[v[0]] = pending
                    pending = None
            elif cur in ("ELEMENT_SOLID", "ELEMENT_SHELL", "ELEMENT_BEAM"):
                v = sl(s, 8)
                if v and len(v) >= 3:
                    elem[cur][v[1]] += 1

print(f"*PART 名称 {len(parts):,} 个\n")
hits = [(pid, t) for pid, t in parts.items() if PAT.search(t)]
print(f"=== 匹配 cart/chondr/menisc/synov/disc/soft/tissue 的部件 ({len(hits)}) ===")
for pid, t in sorted(hits):
    counts = {kw: elem[kw].get(pid, 0) for kw in ("ELEMENT_SOLID", "ELEMENT_SHELL",
                                                 "ELEMENT_BEAM")}
    counts = {k: v for k, v in counts.items() if v}
    print(f"  {pid:>10d}  {counts}  {t}")

print("\n=== 踝/足区段全部部件（pid 8100xxxx–8130xxxx, 8200xxxx–8230xxxx）===")
for pid in sorted(parts):
    if 81000000 <= pid <= 81350000 or 82000000 <= pid <= 82350000:
        counts = {kw: elem[kw].get(pid, 0) for kw in ("ELEMENT_SOLID", "ELEMENT_SHELL",
                                                     "ELEMENT_BEAM")}
        counts = {k: v for k, v in counts.items() if v}
        print(f"  {pid:>10d}  {str(counts):52s} {parts[pid]}")
