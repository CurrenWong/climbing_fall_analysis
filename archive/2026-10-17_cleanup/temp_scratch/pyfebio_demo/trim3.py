"""自动循环精简 MEMORY.md 踝关节段至磁盘 ≤10240。"""
from pathlib import Path

WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
M = WS / "MEMORY.md"
LIMIT = 10240

s = M.read_text(encoding="utf-8")
lines = s.split("\n")
i0 = next(i for i, l in enumerate(lines) if l.startswith("- **接手入口**"))
i1 = next(i for i in range(i0, len(lines)) if lines[i].startswith("- **详见**"))

# 依次尝试的“减法”（每步都保留验证锚点：ANKLE_LIG_TODO / TRUNC_MM=60 / 排序=负结果 / 图 4）
STEPS = [
    ("管线行尾", 2, lambda t: t.replace(
        "）；几何源 THUMS V7.1（E 已是 MPa 零换算）", "）")),
    ("160倍明细", 5, lambda t: t.replace(
        "（32488→9109 单元，**160 倍**）", "（**160 倍**）")),
    ("负结果解释", 6, lambda t: t.replace(
        "——纯压缩下胫骨朝固定足骨**靠近** ⇒ 胫-足韧带**全缩短**",
        "——纯压缩下胫骨朝固定足骨**靠近**")),
    ("韧带清单", 3, lambda t: t.replace(
        "（胫跟=`81101202`、胫舟=`81101203`、", "（胫跟=`81101202`、")),
    ("未完成明细", 7, lambda t: t.replace(
        "subtalar 接触面非 2-流形；**下一步 = 拿吴恺「图 4」原始数值**",
        "**下一步 = 拿吴恺「图 4」原始数值**")),
]

cur = lines[i0:i1 + 1]
for name, idx, fn in STEPS:
    if idx < len(cur):
        cur[idx] = fn(cur[idx])
    out = "\n".join(lines[:i0] + cur + lines[i1 + 1:])
    M.write_text(out, encoding="utf-8")
    sz = M.stat().st_size
    print(f"  去掉 {name:12s} -> disk {sz:,} (余量 {LIMIT - sz:+,})")
    if sz <= LIMIT:
        break

print(f"\n最终 disk = {M.stat().st_size:,} / {LIMIT}")
for i, l in enumerate(cur):
    print(f"  {i}: {l[:100]}")
