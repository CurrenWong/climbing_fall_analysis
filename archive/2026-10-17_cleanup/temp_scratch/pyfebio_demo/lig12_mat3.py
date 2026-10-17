"""L3b: 健壮解析 *MAT 卡（处理 _TITLE 变体的标题行）+ 全部 SECTION_SHELL。"""
import zipfile
from collections import defaultdict
from pathlib import Path

ZIP = Path(r"D:\Project\climbing_fall_analysis\paper\AM50_V71_Occupant.zip")
INNER = "AM50_V71_Occupant/THUMS_model/THUMS_AM50_V71_Occupant_202411.k"

LIG = {81100701: "TALOFIBULARE_POSTERIUS", 81100801: "CALCANEOFIBULARE",
       81100901: "TALOFIBULARE_ANTERIUS", 81101001: "TALOCALCANEUM_POSTERIUS",
       81101101: "TIBIOFIBULARE", 81101201: "MEDIALE1(胫距后深)",
       81101202: "MEDIALE2(胫跟)", 81101203: "MEDIALE3(胫舟)",
       81101204: "MEDIALE4(胫距前深)", 81101401: "LONG_PLANTAR",
       81100100: "COLLAT_FIBULAR_TISSUES", 81100200: "COLLAT_TIBIAL_TISSUES",
       81101410: "CAPITIS_FIBULAE"}


def w10(s):
    return [s[i:i + 10].strip() for i in range(0, len(s) - len(s) % 10, 10)]


mats = defaultdict(list)
secs = defaultdict(list)
cur, buf = None, []
cur_id, cur_lines = None, []


def close_mat():
    global cur_id, cur_lines
    if cur_id is not None and cur_lines:
        mats[cur_id] = ([("KW", cur_lines[0][0])] + 
                        [("", l) for l in cur_lines[1:]]) if False else cur_lines
    cur_id, cur_lines = None, []


def close_sec():
    global cur_id, cur_lines
    if cur_id is not None and cur_lines:
        secs[cur_id] = cur_lines
    cur_id, cur_lines = None, []


with zipfile.ZipFile(ZIP) as z, z.open(INNER) as f:
    for raw in f:
        line = raw.decode("latin-1").rstrip("\r\n")
        st = line.strip()
        if st.startswith("*"):
            if cur and cur.startswith("MAT_"):
                close_mat()
            elif cur == "SECTION_SHELL":
                close_sec()
            kw = st.split()[0].upper().lstrip("*") if st.split() else ""
            cur = kw
            cur_id, cur_lines = None, []
            continue
        if not st or st.startswith("$"):
            continue

        if cur and cur.startswith("MAT_"):
            # 还没定 ID：逐行试第一字段
            if cur_id is None:
                v = w10(st)
                if v and v[0]:
                    try:
                        cur_id = int(v[0])
                        cur_lines.append((cur, st))
                    except ValueError:
                        cur_lines.append(("(title)", st))   # _TITLE 的标题行
            else:
                cur_lines.append(("", st))

        elif cur == "SECTION_SHELL":
            buf.append(st)
            if len(buf) == 2:
                v = w10(buf[0])
                try:
                    secs[int(v[0])] = list(buf)
                except (ValueError, IndexError):
                    pass
                buf = []

close_mat() if cur and cur.startswith("MAT_") else None
close_sec() if cur == "SECTION_SHELL" else None

print(f"解析到 *MAT 卡 {len(mats)} 个，*SECTION_SHELL {len(secs)} 个\n")
print("=" * 74)
for pid, nm in LIG.items():
    print(f"\n### {pid}  {nm}")
    if pid in mats:
        for tag, l in mats[pid]:
            print(f"  MAT {tag:34s} |{l}|")
    else:
        print("  MAT  (无)")
    if pid in secs:
        for l in secs[pid]:
            print(f"  SEC                        |{l}|")
    else:
        print("  SEC  (无 SECTION_SHELL)")

print("\n\n=== 所有 *MAT_FABRIC* 的 id 列表（韧带可能都用这个）===")
fab = sorted(k for k, v in mats.items()
             if v and any("FABRIC" in t for t, _ in v))
print(f"  共 {len(fab)}: {fab[:30]}")
print("\n=== 所有 *MAT_SIMPLIFIED_RUBBER 的 id ===")
rub = sorted(k for k, v in mats.items()
             if v and any("SIMPLIFIED_RUBBER" in t for t, _ in v))
print(f"  共 {len(rub)}: {rub[:30]}")
