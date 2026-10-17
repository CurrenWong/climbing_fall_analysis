"""诊断 v3：全部按 8 字符定宽切片，统计「每个 pid 的单元数」，并用 *PART 标题反查名字。

目的：不再靠猜 pid，直接把 pid -> 名称 -> 单元数 的真相打出来。
"""
import re
import zipfile
from collections import Counter, defaultdict

ZIP = r"D:\Project\climbing_fall_analysis\paper\AM50_V71_Occupant.zip"
INNER = "AM50_V71_Occupant/THUMS_model/THUMS_AM50_V71_Occupant_202411.k"

RE_KW = re.compile(r"^\*([A-Z_0-9]+)")
ELEM_KW = {"ELEMENT_SOLID", "ELEMENT_SHELL", "ELEMENT_SHELL_THICKNESS", "ELEMENT_BEAM"}

parts = {}                  # pid -> title
elem_pids = defaultdict(Counter)     # kw -> Counter(pid)
n_lines = Counter()
samples = defaultdict(list)
n_parse_fail = Counter()
FIELD = 8


def fields(s: str) -> list[int] | None:
    vals = []
    for i in range(0, len(s) - len(s) % FIELD, FIELD):
        f = s[i : i + FIELD].strip()
        if not f:
            continue
        try:
            vals.append(int(f))
        except ValueError:
            return None
    return vals


cur = None
pending_title = None
buf = ""
nbytes = 0

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
                if cur == "PART":
                    pending_title = None
                continue

            if cur == "PART":
                if pending_title is None:
                    pending_title = st
                else:
                    v = fields(st)
                    if v:
                        parts[v[0]] = pending_title
                    pending_title = None
                continue

            if cur in ELEM_KW:
                n_lines[cur] += 1
                v = fields(s)
                if v is None or len(v) < 3:
                    n_parse_fail[cur] += 1
                    if len(samples[cur]) < 5:
                        samples[cur].append(s)
                    continue
                elem_pids[cur][v[1]] += 1

print(f"扫描 {nbytes/1024/1024:.1f} MB   *PART 解析出 {len(parts):,} 个部件名\n")

BONE_HINT = re.compile(
    r"tibia|fibula|femur|patella|talus|calcan|navicular|cuboid|cuneiform|"
    r"metatars|phalang|hallux|^R_|^L_",
    re.I,
)

for kw in ("ELEMENT_SOLID", "ELEMENT_SHELL", "ELEMENT_SHELL_THICKNESS", "ELEMENT_BEAM"):
    print(f"=== *{kw} ===  行 {n_lines[kw]:,}  解析失败 {n_parse_fail[kw]:,}")
    tot = sum(elem_pids[kw].values())
    print(f"  解析出 {tot:,} 个单元，分布在 {len(elem_pids[kw])} 个 pid")
    # 只打与骨/韧带相关的
    hits = []
    for pid, cnt in elem_pids[kw].most_common():
        title = parts.get(pid, "<无*PART标题>")
        if BONE_HINT.search(title) or "LIG" in title.upper() or "BONE" in title.upper():
            hits.append((pid, cnt, title))
    print(f"  骨/韧带相关 pid: {len(hits)} 个")
    for pid, cnt, title in hits[:60]:
        print(f"    {pid:>10d}  {cnt:>7,}  {title}")
    if len(hits) > 60:
        print(f"    ... 其余 {len(hits)-60} 个略")
    if samples[kw]:
        print("  解析失败样本:")
        for s in samples[kw]:
            print(f"    len={len(s)} {s!r}")
    print()
