"""THUMS .k 部件发现：正确解析后列出 pid -> *PART 名称 -> 单元数。

实测格式（务必记住，两种宽度并存）：
  *NODE               : I8 + 3*F16   （有空格，空白分隔即可）
  *ELEMENT_SOLID      : **8 字符定宽，无空格**  eid,pid,n1..n8   （行=80 字符）
  *ELEMENT_SHELL      : **8 字符定宽，无空格**  eid,pid,n1..n4   （行=48 字符）
  *ELEMENT_BEAM       : **8 字符定宽**          eid,pid,n1,n2,...
  *ELEMENT_SHELL_THICKNESS : 厚度值（非单元），10 字符浮点
  *PART 数值卡         : **10 字符定宽**  pid,secid,mid,...      （行=90 字符）
"""
import re
import sys
import zipfile
from collections import Counter, defaultdict

ZIP = r"D:\Project\climbing_fall_analysis\paper\AM50_V71_Occupant.zip"
INNER = "AM50_V71_Occupant/THUMS_model/THUMS_AM50_V71_Occupant_202411.k"

RE_KW = re.compile(r"^\*([A-Z_0-9]+)")
ELEM_W8 = {"ELEMENT_SOLID", "ELEMENT_SHELL", "ELEMENT_BEAM"}


def slice_fields(s: str, width: int) -> list[int] | None:
    vals = []
    for i in range(0, len(s) - len(s) % width, width):
        f = s[i : i + width].strip()
        if not f:
            continue
        try:
            vals.append(int(f))
        except ValueError:
            return None
    return vals


def main() -> int:
    parts: dict[int, str] = {}
    elem_pids: dict[str, Counter] = defaultdict(Counter)
    n_lines = Counter()
    fails = Counter()
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
                        pending_title = st          # 第 1 行 = 部件标题
                    else:
                        v = slice_fields(st, 10)    # ★ 10 字符
                        if v:
                            parts[v[0]] = pending_title
                        else:
                            fails["PARTカード"] += 1
                        pending_title = None
                    continue

                if cur in ELEM_W8:
                    n_lines[cur] += 1
                    v = slice_fields(s, 8)          # ★ 8 字符
                    if v is None or len(v) < 3:
                        fails[cur] += 1
                        continue
                    elem_pids[cur][v[1]] += 1

    print(f"扫描 {nbytes/1024/1024:.1f} MB")
    print(f"*PART: 解析出 {len(parts):,} 个部件名，失败 {fails['PARTカード']}")
    for kw in sorted(ELEM_W8):
        print(f"*{kw}: 行 {n_lines[kw]:,}  解析失败 {fails[kw]:,}  "
              f"pid 数 {len(elem_pids[kw])}")
    print()

    # ---- 目标名称模式 ----
    PAT = re.compile(
        r"(tibia|fibula|femur|patella|talus|calcan|navicular|cuboid|cuneiform|"
        r"metatars|phalang|lig_|ligament|foot)", re.I)

    rows = []
    for kw in sorted(ELEM_W8):
        for pid, cnt in elem_pids[kw].items():
            title = parts.get(pid, "<未命名>")
            if PAT.search(title):
                rows.append((title.lower(), pid, kw, cnt, title))

    rows.sort()
    print("=== 骨 / 韧带 / 足部 相关部件（按名称排序）===")
    print(f"{'pid':>10s}  {'keyword':24s} {'elems':>8s}  *PART 名称")
    for _, pid, kw, cnt, title in rows:
        print(f"{pid:>10d}  {kw:24s} {cnt:>8,}  {title}")
    print(f"\n共 {len(rows)} 条")

    print("\n=== 全部未命名 pid（前 20）===")
    named = set(parts)
    unn = [(kw, pid, c) for kw in ELEM_W8 for pid, c in elem_pids[kw].items()
           if pid not in named]
    print(f"  共 {len(unn)} 个")
    for kw, pid, c in unn[:20]:
        print(f"  {pid:>10d}  {kw:24s} {c:>8,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
