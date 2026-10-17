"""dump *PART 块与 *SECTION_SHELL / *MAT 原始行（看清真实格式）。"""
import zipfile
from pathlib import Path

ZIP = Path(r"D:\Project\climbing_fall_analysis\paper\AM50_V71_Occupant.zip")
INNER = "AM50_V71_Occupant/THUMS_model/THUMS_AM50_V71_Occupant_202411.k"

cur = None
buf = []
shown = 0
with zipfile.ZipFile(ZIP) as z, z.open(INNER) as f:
    for raw in f:
        line = raw.decode("latin-1").rstrip("\r\n")
        st = line.strip()
        if st.startswith("*"):
            kw = st.split()[0].upper().lstrip("*")
            # 打印刚结束的 PART 块
            if cur == "PART" and any("811" in l for l in buf) and shown < 3:
                print(f"--- *PART block (含 811) ---")
                for l in buf[:8]:
                    print(f"    {len(l):3d} |{l}|")
                shown += 1
            cur = kw
            buf = []
            if kw == "PART":
                buf.append("<<*PART>>")
            continue
        buf.append(line)
        if len(buf) > 12:
            buf.pop(0)

print("\n=== *SECTION_SHELL 块（前 3 个）===")
cur, buf, shown = None, [], 0
with zipfile.ZipFile(ZIP) as z, z.open(INNER) as f:
    for raw in f:
        line = raw.decode("latin-1").rstrip("\r\n")
        st = line.strip()
        if st.startswith("*"):
            kw = st.split()[0].upper().lstrip("*")
            if cur == "SECTION_SHELL" and shown < 3:
                print(f"--- SECTION_SHELL ---")
                for l in buf[:6]:
                    print(f"    {len(l):3d} |{l}|")
                shown += 1
            cur = kw
            buf = []
            continue
        if cur == "SECTION_SHELL":
            buf.append(line)
        if len(buf) > 8:
            buf.pop(0)
