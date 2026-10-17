"""探测缺失部件的单元类型（solid/shell），为扩充 TARGETS 做准备。"""
import re
import zipfile
from collections import defaultdict
from pathlib import Path

ZIP = Path(r"D:\Project\climbing_fall_analysis\paper\AM50_V71_Occupant.zip")
INNER = "AM50_V71_Occupant/THUMS_model/THUMS_AM50_V71_Occupant_202411.k"

WANT = set()
for p in (81001400, 81001401, 81001500, 81001501, 81001600, 81001601,
          81001700, 81001701, 81001800, 81001801,
          81001900, 81001901, 81002000, 81002001,
          81101301, 81101302, 81101303, 81101304, 81101305, 81101306,
          81101307, 81101308, 81101309,
          81101408, 81101409, 81101410):
    WANT.add(p)

ELEM = {"ELEMENT_SOLID": 8, "ELEMENT_SHELL": 4, "ELEMENT_BEAM": 2}
cnt = defaultdict(lambda: defaultdict(int))
lastkw = None
with zipfile.ZipFile(ZIP) as z, z.open(INNER) as f:
    for raw in f:
        line = raw.decode("latin-1").rstrip("\r\n")
        if line.startswith("*"):
            kw = line.split()[0].upper().lstrip("*") if line.split() else ""
            lastkw = kw
            # 单元块带 pid 选项也可能在 header 上；此处用定宽切片取 pid
            continue
        if lastkw in ELEM:
            w = 8
            s = line
            if len(s) < w * 2:
                continue
            try:
                eid = int(s[0:w])
                pid = int(s[w:w * 2])
            except ValueError:
                continue
            if pid in WANT:
                cnt[pid][lastkw] += 1

print("pid        类型           单元数")
for pid in sorted(cnt):
    for k, n in cnt[pid].items():
        print(f"  {pid}  {k:16s} {n:6d}")
print(f"\n命中 {len(cnt)}/{len(WANT)} 个部件")
missing = WANT - set(cnt)
if missing:
    print("未命中:", sorted(missing))
