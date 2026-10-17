"""诊断：THUMS .k 里 *ELEMENT_SOLID / *ELEMENT_SHELL 数据行的真实格式。

按块统计「每行 token 数」直方图，并抓中段样本。
若直方图出现非预期的 token 数，说明字段粘连（定宽），需改解析策略。
"""
import re
import zipfile
from collections import Counter, defaultdict

ZIP = r"D:\Project\climbing_fall_analysis\paper\AM50_V71_Occupant.zip"
INNER = "AM50_V71_Occupant/THUMS_model/THUMS_AM50_V71_Occupant_202411.k"

RE_KW = re.compile(r"^\*([A-Z_0-9]+)")
INTEREST = {"ELEMENT_SOLID", "ELEMENT_SHELL"}

tok_hist = defaultdict(Counter)      # kw -> Counter(token count)
len_hist = defaultdict(Counter)      # kw -> Counter(line length)
samples = defaultdict(list)          # kw -> mid-block samples
nlines = Counter()
pids_seen = defaultdict(Counter)     # kw -> Counter(pid)
pid_first_title = {}

cur = None
buf = ""
nbytes = 0
n_data = Counter()

with zipfile.ZipFile(ZIP) as z, z.open(INNER) as f:
    while True:
        chunk = f.read(8 * 1024 * 1024)
        if not chunk:
            break
        nbytes += len(chunk)
        buf += chunk.decode("latin-1")
        lines = buf.split("\n")
        buf = lines.pop()
        for raw in lines:
            s = raw.rstrip("\r")
            st = s.strip()
            if not st or st[0] == "$":
                continue
            if st[0] == "*":
                m = RE_KW.match(st.upper())
                cur = m.group(1) if m else None
                continue
            nlines[cur] += 1
            if cur in INTEREST:
                t = st.split()
                tok_hist[cur][len(t)] += 1
                len_hist[cur][len(s)] += 1
                n_data[cur] += 1
                if len(t) >= 2:
                    try:
                        pids_seen[cur][int(t[1])] += 1
                    except ValueError:
                        pass
                # 收集第 2000~2010 行做样本
                if 2000 <= n_data[cur] <= 2010:
                    samples[cur].append(s)

print(f"总扫描 {nbytes/1024/1024:.1f} MB\n")
for kw in INTEREST:
    print(f"=== *{kw} ===")
    print(f"  数据行 {n_data[kw]:,}")
    print(f"  token 数直方图: {dict(sorted(tok_hist[kw].items()))}")
    print(f"  行长度 top6: {sorted(len_hist[kw].items(), key=lambda x: -x[1])[:6]}")
    top = pids_seen[kw].most_common(8)
    print(f"  pid top8: {[(p, c) for p, c in top]}")
    # 目标 pid 是否出现
    for want in (81000900, 81000700, 81000800, 81000801, 81000600, 81000601,
                 81001100, 81001300):
        print(f"    pid {want}: {pids_seen[kw].get(want, 0):,} 行")
    print("  中段样本:")
    for s in samples[kw][:6]:
        print(f"    len={len(s):3d}  {s!r}")
    print()

print("=== 其它可能含单元的关键字（数据行数）===")
for kw, c in nlines.most_common(30):
    if kw and ("ELEMENT" in kw or "SET" in kw or "PART" in kw or "NODE" in kw):
        print(f"  {kw:38s} {c:>10,}")
