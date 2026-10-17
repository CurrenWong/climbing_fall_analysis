"""查当前 .feb 里的韧带杆：数量 + ATFL 配对是否已分布（验证修复是否生效）。"""
import io
import re
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
feb = ROOT / "temp/thums/febio_lig/ankle_lig_WK.feb"
t = feb.read_text(encoding="latin-1", errors="replace")

# 韧带单元行（line2 元素）
m = re.search(r'<Elements[^>]*type="line2"[^>]*>(.*?)</Elements>', t, re.S)
if not m:
    print("未找到 line2 元素段；列出所有 Elements 段的 type 与行数：")
    for mm in re.finditer(r'<Elements([^>]*)>(.*?)</Elements>', t, re.S):
        attrs = mm.group(1).strip()[:90]
        n = len([x for x in mm.group(2).split("\n") if x.strip()])
        print(f"  {attrs:<92} {n} 行")
else:
    body = m.group(1)
    rows = [[int(x) for x in l.split(",")] for l in body.split("\n") if l.strip()]
    print(f".feb 中 line2（韧带）单元 {len(rows)} 个")
    # 与 ligaments.npz 对比
    z = np.load(ROOT / "temp/thums/ligaments.npz", allow_pickle=True)
    names = [str(x) for x in z["names"]]
    print(f"ligaments.npz: 杆总数 {int(np.sum(z['n_truss']))}，各韧带 "
          f"{[int(x) for x in z['n_truss']]}")
    print(f"  names = {names}")
    # 看前 8 个单元
    print("前 8 个线单元（node1,node2,...）：")
    for r in rows[:8]:
        print("   ", r[:8])
    print(f"\n★ 若杆总数 = {int(np.sum(z['n_truss']))} ⇒ .feb 用的就是新 ligaments.npz")

# ATFL 的 truss（从 npz 直接看）
z = np.load(ROOT / "temp/thums/ligaments.npz", allow_pickle=True)
i = [k for k, n in enumerate(z["names"]) if "ATFL" in str(n)]
if i:
    k = i[0]
    n = int(z["n_truss"][k])
    tr = z["truss"][k][:n]
    print(f"\nATFL（npz）n_truss={n}:")
    for a, b in tr:
        print(f"    {int(a)} <-> {int(b)}")
    print(f"  不同 B 端点数 = {len({int(b) for _, b in tr})} "
          f"（=1 则是退化扇形）")
