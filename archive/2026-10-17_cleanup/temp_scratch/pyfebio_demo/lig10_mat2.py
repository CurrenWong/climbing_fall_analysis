"""L3: 正确解析 *PART(pid,secid,mid) -> *SECTION_SHELL(厚度) + *MAT(材料)。"""
import zipfile
from collections import defaultdict
from pathlib import Path

ZIP = Path(r"D:\Project\climbing_fall_analysis\paper\AM50_V71_Occupant.zip")
INNER = "AM50_V71_Occupant/THUMS_model/THUMS_AM50_V71_Occupant_202411.k"

WANT = {81100100, 81100200, 81100701, 81100801, 81100901, 81101001,
        81101101, 81101201, 81101202, 81101203, 81101204, 81101401,
        81101410, 81101301}

pid2 = {}                      # pid -> (secid, mid)
secs = defaultdict(list)       # secid -> [lines]
mats = defaultdict(list)       # mid -> [(kw, line)]

cur, buf = None, []
curmat_id = None


def w10(s):
    out = []
    for i in range(0, len(s) - len(s) % 10, 10):
        f = s[i:i + 10].strip()
        out.append(f)
    return out


with zipfile.ZipFile(ZIP) as z, z.open(INNER) as f:
    for raw in f:
        line = raw.decode("latin-1").rstrip("\r\n")
        st = line.strip()
        if st.startswith("*"):
            kw = st.split()[0].upper().lstrip("*") if st.split() else ""
            cur, buf = kw, []
            curmat_id = None
            continue
        if not st or st.startswith("$"):
            continue

        if cur == "PART":
            buf.append(st)
            if len(buf) == 2:                       # [title, 数值卡]
                v = w10(buf[1])
                try:
                    pid, secid, mid = int(v[0]), int(v[1]), int(v[2])
                except (ValueError, IndexError):
                    continue
                if pid in WANT:
                    pid2[pid] = (secid, mid)
                buf = []

        elif cur == "SECTION_SHELL":
            buf.append(st)
            if len(buf) == 2:                       # [secid行, 厚度行]
                v = w10(buf[0])
                try:
                    secs[int(v[0])] = [buf[0], buf[1]]
                except (ValueError, IndexError):
                    pass
                buf = []

        elif cur.startswith("MAT_"):
            if curmat_id is None:
                try:
                    curmat_id = int(w10(st)[0])
                except (ValueError, IndexError):
                    curmat_id = -1
                mats[curmat_id].append((cur, st))
            else:
                mats[curmat_id].append(("", st))

print("=== pid -> (secid, mid) ===")
for p in sorted(pid2):
    print(f"  {p} -> {pid2[p]}")

seen = set()
for p in sorted(pid2):
    secid, mid = pid2[p]
    if mid in seen:
        continue
    seen.add(mid)
    print(f"\n{'='*70}\n  pid={p}  secid={secid}  mid={mid}")
    print("  --- *MAT ---")
    for kw, l in mats.get(mid, [("", "(无 MAT 卡)")]):
        print(f"    {kw:34s} |{l}|")
    print("  --- *SECTION_SHELL ---")
    for l in secs.get(secid, ["(无 SECTION_SHELL)"]):
        print(f"    |{l}|")
