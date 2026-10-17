"""MEMORY.md 踝关节段等量替换：加入「刚度 2 倍不确定度」这条关键警告。"""
from pathlib import Path

WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
M = WS / "MEMORY.md"
LIMIT = 10240
L = M.read_text(encoding="utf-8").split("\n")
i0 = next(i for i, l in enumerate(L) if l.startswith("- **接手入口**"))
i1 = next(i for i in range(i0, len(L)) if L[i].startswith("- **详见**"))

NEW = [
    "- **接手入口**：`FE_PIPELINE_TODO.md`「🧭 接手入口」节；**刚度问题**→`CONV_TODO.md`",
    "- ⚠️ **刚度有 2 倍量级不确定度（未归因）**：同几何/材料/边界下 tet4 **259.3 N** vs",
    "  penta6 **116.9 N**（DISP=0.1472）。五假设已排除（单元格式/锁死/几何/剔除/接触面积），",
    "  最尖锐候选=**接触面定向**。⇒ 600 N 档**应力 2.96 MPa / 面积 287 mm² 双靶可用**，",
    "  **载荷-位移关系不可信**",
    "- **管线**：C/B/A 全 ✅（单关节 `rc=0`/16.6 s，436.8 N，软骨 p99 2.58/4.27",
    "  → 外推 600 N ≈3.4 MPa，对吴恺 3.74）；源 THUMS V7.1",
    "- **韧带轮**（→`ANKLE_LIG_TODO.md`）：4 条已定位（胫跟=`81101202`、跟腓=`81100801`、",
    "  下胫腓=`81101101`；EA 510/510/**75.4**、厚 2.0）；`TRUNC_MM=60` 截断胫腓骨",
    "  （**160 倍**，`rc=0`/5.4 s）",
    "- ❌ **排序=负结果**：ΔL 0.001~0.08 vs 吴恺 1.92~5.24 mm；胫跟名次颠倒",
    "- **详见**：`memory/2026-10-03/ankle-fe-modeling.md` §六–§十",
]

out = "\n".join(L[:i0] + NEW + L[i1 + 1:])
M.write_text(out, encoding="utf-8")
sz = M.stat().st_size
print(f"disk = {sz:,} / {LIMIT}  (余量 {LIMIT - sz:+,})")
if sz > LIMIT:
    print("!! 超限，需再裁")

# 复核
L2 = M.read_text(encoding="utf-8").split("\n")
print("\n--- 新段 ---")
for i in range(i0, i0 + len(NEW)):
    print(f"{i+1:4d}| {L2[i][:96]}")
