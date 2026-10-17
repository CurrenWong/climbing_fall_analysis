"""C1: 抽出并查看 THUMS 骨材料卡的内容。
"""
import zipfile
from pathlib import Path

ZIP = Path(r"D:\Project\climbing_fall_analysis\paper\AM50_V71_Occupant.zip")
OUT = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
OUT.mkdir(parents=True, exist_ok=True)

NAMES = [
    "AM50_V71_Occupant/THUMS_model/mat_bone_AM50_V71_Occ_fracture.k",
    "AM50_V71_Occupant/THUMS_model/mat_bone_AM50_V71_Occ_no_fracture.k",
]

with zipfile.ZipFile(ZIP) as z:
    for n in NAMES:
        try:
            data = z.read(n)
        except KeyError:
            print(f"缺失 {n}")
            continue
        dst = OUT / Path(n).name
        dst.write_bytes(data)
        print(f"\n{'='*80}\n{n}  ({len(data):,} bytes) -> {dst.name}\n{'='*80}")
        txt = data.decode("latin-1")
        lines = txt.split("\n")
        print(f"总行数 {len(lines)}")
        # 打印所有关键字行
        kws = [(i, l.strip()) for i, l in enumerate(lines) if l.strip().startswith("*")]
        print(f"\n--- 关键字（{len(kws)} 个）---")
        for i, l in kws:
            print(f"  L{i:<6d} {l}")
        print("\n--- 前 90 行原样 ---")
        for i, l in enumerate(lines[:90]):
            print(f"  {i:4d}| {l.rstrip()}")
