"""MEMORY 最后一次压缩（紧凑化今天新增的三行踝段）。"""
import io

MEM = r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md"
m = io.open(MEM, encoding="utf-8").read()
n0 = len(m.encode("utf-8"))
REPS = [
    ("（`[ERR-20261006-003]`；三证：`.feb` 写 z=−DISP ｜ 顶面 dz 朝距骨 ｜ 接触 149.6 vs 0 N）。\n"
     "  **上轮 `[ERR-20261005-009]` 判反、已作废**；`LOAD_SIGN` 对 `dispx` 无效。",
     "（`[ERR-20261006-003]`；三证：`.feb` 写 z=−DISP ｜ dz 朝距骨 ｜ 接触 149.6 vs 0 N）。\n"
     "  **上轮 `[ERR-20261005-009]` 判反、已作废**；`LOAD_SIGN` 对 `dispx` 无效。"),
    ("- ✅ **真·压缩曲线**（`[LRN-079]`，全 NORMAL）：+0.25/+0.5/+1.0/**+1.6**/+2.5 mm →\n"
     "  124.0/230.7/418.1/**592.6**/615.4 N ⇒ **600 N 在 `DISP_MM≈+1.9 mm`**\n"
     "  ⚠️ `[LRN-072]` 的 −132…−1017 N 是**拉伸**、**已作废**",
     "- ✅ **真·压缩曲线**（`[LRN-079]`）：+0.25/0.5/1.0/**1.6**/2.5 mm → 124.0/230.7/418.1/**592.6**/615.4 N\n"
     "  ⇒ **600 N 在 `DISP_MM≈+1.9 mm`**；⚠️ `[LRN-072]` 的 −132…−1017 N 是**拉伸**、已作废"),
    ("- 🧭 **载荷路径**（`[LRN-080]`）：压缩下**韧带 77% / 接触 28%**（距骨全固定对照）；\n"
     "  `FIX_TALUS=transverse` 时几乎全走韧带。**路径对照必须 `FIX_TALUS=all`**（否则近奇异）。\n"
     "  `shell` 模式用固定 `E=260`、不用 `a_cross`",
     "- 🧭 **载荷路径**（`[LRN-080]`）：压缩下**韧带 77% / 接触 28%**；**对照必须 `FIX_TALUS=all`**"
     "（否则近奇异）。`shell` 用固定 `E=260`、不用 `a_cross`"),
]
for old, new in REPS:
    if old in m:
        m = m.replace(old, new, 1)
    else:
        print(f"  ⚠ 未命中：{old[:44]}")
io.open(MEM, "w", encoding="utf-8").write(m)
n1 = len(m.encode("utf-8"))
print(f"MEMORY.md {n0} -> {n1} B  {'✅ ≤10240' if n1 <= 10240 else f'❌ 超 {n1-10240}'}")
