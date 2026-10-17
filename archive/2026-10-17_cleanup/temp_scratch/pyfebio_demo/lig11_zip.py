"""查 zip 内全部文件 + 主 .k 的 *INCLUDE 指令。"""
import re
import zipfile
from pathlib import Path

ZIP = Path(r"D:\Project\climbing_fall_analysis\paper\AM50_V71_Occupant.zip")
with zipfile.ZipFile(ZIP) as z:
    infos = z.infolist()
    print(f"zip 内 {len(infos)} 个条目\n")
    for i in infos:
        if not i.is_dir():
            print(f"  {i.file_size:>12,}  {i.filename}")
    # 主 .k 的 INCLUDE
    inner = [i.filename for i in infos if i.filename.endswith(".k")]
    print(f"\n.k 文件 {len(inner)} 个:")
    for n in inner:
        print("   ", n)

main = "AM50_V71_Occupant/THUMS_model/THUMS_AM50_V71_Occupant_202411.k"
with zipfile.ZipFile(ZIP) as z, z.open(main) as f:
    incs = []
    for raw in f:
        line = raw.decode("latin-1").rstrip("\r\n")
        if line.strip().upper().startswith("*INCLUDE"):
            incs.append(line.strip())
    print(f"\n主 .k 内 *INCLUDE {len(incs)} 条:")
    for l in incs[:30]:
        print("   ", l)
