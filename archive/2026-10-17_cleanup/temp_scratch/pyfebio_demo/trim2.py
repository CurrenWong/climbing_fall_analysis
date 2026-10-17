"""精简 MEMORY.md 踝关节段（磁盘口径 ≤10240），保留验证锚点。"""
from pathlib import Path

WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
M = WS / "MEMORY.md"

s = M.read_text(encoding="utf-8")
lines = s.split("\n")
i0 = next(i for i, l in enumerate(lines) if l.startswith("- **接手入口**"))
i1 = next(i for i in range(i0, len(lines)) if lines[i].startswith("- **详见**"))

cur = lines[i0:i1 + 1]
print("当前段:")
for i, l in enumerate(cur):
    print(f"  {i}: {l}")

# 精简（保留 ANKLE_LIG_TODO / TRUNC_MM=60 / 排序=负结果 / 图 4 这些锚点）
REPL = {
    8: "- **详见**：`memory/2026-10-03/ankle-fe-modeling.md` §六–§十",
    7: "- **未完成**：subtalar 接触面非 2-流形；**下一步 = 拿吴恺「图 4」原始数值**",
    5: "  **跑通 `rc=0`/5.4 s**，关键 **`TRUNC_MM=60` 截断胫腓骨**（32488→9109，**160 倍**）",
    6: "- ❌ **排序=负结果**：ΔL **0.001~0.08** vs 吴恺 1.92~5.24 mm；**胫跟名次颠倒**——纯压缩下胫骨朝固定足骨**靠近** ⇒ 胫-足韧带**全缩短**",
    2: "  → 外推 600 N 胫侧 ≈3.4 MPa，**与吴恺 3.74 MPa 吻合**）；源 THUMS V7.1",
}
for k, v in REPL.items():
    if k < len(cur):
        cur[k] = v

out = "\n".join(lines[:i0] + cur + lines[i1 + 1:])
M.write_text(out, encoding="utf-8")
print(f"\ndisk = {M.stat().st_size:,} / 10240  -> "
      f"{'OK' if M.stat().st_size <= 10240 else 'STILL OVER'}")
print(f"headroom = {10240 - M.stat().st_size:+,}")
