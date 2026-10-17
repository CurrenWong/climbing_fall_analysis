"""按**磁盘口径**（含 CRLF）收敛 MEMORY.md 到 10240 以内。"""
from pathlib import Path

WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
M = WS / "MEMORY.md"

raw = M.read_bytes()
print(f"磁盘 {len(raw):,} bytes | CRLF 数 {raw.count(b'\\r\\n')} | "
      f"LF 数 {raw.count(b'\\n')}")

s = M.read_text(encoding="utf-8")
lines = s.split("\n")

i0 = next(i for i, l in enumerate(lines) if l.startswith("- **接手入口**"))
i1 = next(i for i in range(i0, len(lines)) if lines[i].startswith("- **详见**"))

# 逐步精简（按磁盘口径判断）
FALLBACK = [
    (8, "- **详见**：`memory/2026-10-03/ankle-fe-modeling.md` §六–§十",
        "- **详见**：`memory/2026-10-03/ankle-fe-modeling.md` §六–§十"),
    (7, "- **未完成**：subtalar 接触面非 2-流形；**下一步 = 拿吴恺「图 4」原始数值**（文本抽取不含）",
        "- **未完成**：subtalar 接触面非 2-流形；**下一步 = 拿吴恺「图 4」原始数值**"),
    (1, "  → 外推 600 N 胫侧 ≈3.4 MPa，**与吴恺 3.74 MPa 吻合**）；几何源 THUMS V7.1（E 已是 MPa 零换算）",
        "  → 外推 600 N 胫侧 ≈3.4 MPa，**与吴恺 3.74 MPa 吻合**）；源 THUMS V7.1"),
]

cur = lines[i0:i1 + 1]
for idx, old, new in FALLBACK:
    if idx < len(cur) and cur[idx] == old:
        cur[idx] = new
    out = "\n".join(lines[:i0] + cur + lines[i1 + 1:])
    disk = len(out.replace("\r\n", "\n").replace("\n", "\r\n").encode("utf-8"))
    print(f"  try -> 磁盘 {disk:,}")
    if disk <= 10240:
        break

# 保留原行尾风格写回
txt = "\n".join(lines[:i0] + cur + lines[i1 + 1:])
M.write_bytes(txt.replace("\r\n", "\n").replace("\n", "\r\n").encode("utf-8"))
raw2 = M.read_bytes()
print(f"最终磁盘 {len(raw2):,} / 10240  -> "
      f"{'OK' if len(raw2) <= 10240 else 'STILL OVER'}")
