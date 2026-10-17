"""L3: 韧带的 *PART(mid) -> *MAT 卡 + 壳厚来源（*SECTION_SHELL / SHELL_THICKNESS）。"""
import re
import zipfile
from collections import defaultdict
from pathlib import Path

ZIP = Path(r"D:\Project\climbing_fall_analysis\paper\AM50_V71_Occupant.zip")
INNER = "AM50_V71_Occupant/THUMS_model/THUMS_AM50_V71_Occupant_202411.k"

LIG_PIDS = [81100100, 81100200, 81100701, 81100801, 81100901, 81101001,
            81101101, 81101201, 81101202, 81101203, 81101204, 81101401,
            81101410]
WANT = set(LIG_PIDS)

pid_to_mid = {}          # pid -> (secid, mid)
mat_cards = defaultdict(list)   # mid -> lines
sec_cards = defaultdict(list)   # secid -> lines
thick_by_pid = defaultdict(list)

cur = None
part_buf = []
mat_id = None
sec_id = None

with zipfile.ZipFile(ZIP) as z, z.open(INNER) as f:
    for raw in f:
        line = raw.decode("latin-1").rstrip("\r\n")
        st = line.strip()
        if st.startswith("*"):
            # 收尾上一块
            if cur == "PART" and part_buf:
                v = []
                s = part_buf[0]
                for i in range(0, len(s) - len(s) % 10, 10):
                    fld = s[i:i + 10].strip()
                    if fld:
                        try:
                            v.append(int(fld))
                        except ValueError:
                            v = []
                            break
                if len(v) >= 3 and v[0] in WANT:
                    pid_to_mid[v[0]] = (v[1], v[2])
                part_buf = []
            kw = st.split()[0].upper().lstrip("*") if st.split() else ""
            cur = kw
            if kw.startswith("MAT_"):
                # *MAT_024 / *MAT_ELASTIC 等，第一行是 id
                mat_id = None
            elif kw.startswith("SECTION_SHELL"):
                sec_id = None
            continue
        if not st or st.startswith("$"):
            continue

        if cur == "PART":
            part_buf.append(st)

        elif cur.startswith("MAT_"):
            if mat_id is None:
                # 第一行第一个字段是 mid
                fld = st[:10].strip()
                try:
                    mat_id = int(fld)
                except ValueError:
                    continue
                mat_cards[mat_id].append(("*" + cur, st))
            else:
                mat_cards[mat_id].append(("", st))

        elif cur.startswith("SECTION_SHELL"):
            if sec_id is None:
                fld = st[:10].strip()
                try:
                    sec_id = int(fld)
                except ValueError:
                    continue
                sec_cards[sec_id].append(("*" + cur, st))
            else:
                sec_cards[sec_id].append(("", st))

print("=== pid -> (secid, mid) ===")
for pid in LIG_PIDS:
    print(f"  {pid}: {pid_to_mid.get(pid, 'NOT FOUND')}")

print("\n=== 涉及的 *MAT 卡 ===")
seen = set()
for pid in LIG_PIDS:
    if pid not in pid_to_mid:
        continue
    secid, mid = pid_to_mid[pid]
    if mid in seen:
        continue
    seen.add(mid)
    print(f"\n  --- mid={mid} (用于 pid={pid}) ---")
    for tag, l in mat_cards.get(mid, [("", "(该 mid 无 MAT 卡，可能是 NUL/空材料)")]):
        print(f"    {tag:22s} {l}")
    print(f"    [section secid={secid}]")
    for tag, l in sec_cards.get(secid, [("", "(无 SECTION_SHELL)")]):
        print(f"      {tag:22s} {l}")

print(f"\n\n=== 全部含 LIG / MAT 的卡片名（供参考）===")
allmat = sorted({t for v in mat_cards.values() for t, _ in v if t})
print("  ", allmat[:40])
print(f"  MAT 卡总数 {len(mat_cards)}  |  SECTION_SHELL {len(sec_cards)}")
