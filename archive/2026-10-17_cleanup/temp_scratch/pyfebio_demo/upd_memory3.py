"""精简 + 自动收敛到 10240 以内（按行依次砍冗余）。"""
from pathlib import Path

WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
M = WS / "MEMORY.md"

NEW = [
    "- **接手入口**：`FE_PIPELINE_TODO.md`「🧭 接手入口」节（自包含）；历史 `docs/踝关节FE建模交接.md`",
    "- **管线状态**：**C ✅ / B ✅ / A ✅**（`rc=0`/16.6 s，支反力 **436.8 N**，软骨 p99 **2.58/4.27 MPa**",
    "  → 外推 600 N 胫侧 ≈3.4 MPa，**与吴恺 3.74 MPa 吻合**）；几何源 THUMS V7.1（E 已是 MPa 零换算）",
    "- **韧带轮**（细节→`ANKLE_LIG_TODO.md`）：吴恺 4 条已定位（胫跟=`81101202`、胫舟=`81101203`、",
    "  跟腓=`81100801`、下胫腓=`81101101`；EA 510/510/510/**75.4**、厚 2.0 mm）；",
    "  **跑通 `rc=0`/5.4 s**，关键 **`TRUNC_MM=60` 截断胫腓骨**（32488→9109 单元，**160 倍**）",
    "- ❌ **排序=负结果**：ΔL **0.001~0.08** vs 吴恺 1.92~5.24 mm；**胫跟名次颠倒**——纯压缩下胫骨朝固定足骨**靠近** ⇒ 胫-足韧带**全缩短**",
    "- **未完成**：subtalar 接触面非 2-流形；**下一步 = 拿吴恺「图 4」原始数值**（文本抽取不含）",
    "- **详见**：`memory/2026-10-03/ankle-fe-modeling.md` §六–§十、`2026-10-04/…handover.md` §四",
]

# 按需再砍的备选（依次启用）
FALLBACK = [
    (3, "- **韧带轮**（细节→`ANKLE_LIG_TODO.md`）：吴恺 4 条已定位（胫跟=`81101202`、胫舟=`81101203`、",
        "- **韧带轮**（→`ANKLE_LIG_TODO.md`）：吴恺 4 条已定位（胫跟=`81101202`、胫舟=`81101203`、"),
    (8, "- **详见**：`memory/2026-10-03/ankle-fe-modeling.md` §六–§十、`2026-10-04/…handover.md` §四",
        "- **详见**：`memory/2026-10-03/ankle-fe-modeling.md` §六–§十"),
    (1, "- **管线状态**：**C ✅ / B ✅ / A ✅**（`rc=0`/16.6 s，支反力 **436.8 N**，软骨 p99 **2.58/4.27 MPa**",
        "- **管线**：**C/B/A 全 ✅**（`rc=0`/16.6 s，支反力 **436.8 N**，软骨 p99 **2.58/4.27 MPa**"),
]

s = M.read_text(encoding="utf-8")
lines = s.split("\n")
i0 = next(i for i, l in enumerate(lines) if l.startswith("- **交接文档**"))
i1 = next(i for i in range(i0, len(lines)) if lines[i].startswith("- **详见**"))
head, tail = lines[:i0], lines[i1 + 1:]

cur = list(NEW)
for idx, old, new in FALLBACK:
    out = "\n".join(head + cur + tail)
    if len(out.encode("utf-8")) <= 10240:
        break
    cur[idx] = new
out = "\n".join(head + cur + tail)
b = len(out.encode("utf-8"))
print(f"{b:,} / 10240  (余量 {10240 - b:+,})")
if b <= 10240:
    M.write_text(out, encoding="utf-8")
    print("已写入")
else:
    print("仍超限，需砍别处")
