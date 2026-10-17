"""MEMORY 削 ~165 B：把两条已在 reference 文件里的小程序 gotcha 从主表移除（先验证存在）。"""
from pathlib import Path
import re

WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
ROOT = Path(r"D:\Project\climbing_fall_analysis")
M = WS / "MEMORY.md"
REF = WS / "references/seclog-miniprogram-gotchas.md"

ref = REF.read_text(encoding="utf-8") if REF.exists() else ""
print(f"reference 文件: {REF.exists()}  {len(ref):,} chars")
for key in ("88rpx", "babel", "min-height"):
    print(f"  [{'OK' if key in ref else 'XX'}] reference 含 {key}")

s = M.read_text(encoding="utf-8")
CAND = [
    "| 小程序 input 比 picker 矮 | 原生组件不响应 `line-height` ⇒ `min-height: 88rpx` |\n",
    "| miniprogram-ci `Unexpected token .` | babel 不支持 `?.` ⇒ `(a && a.b)` |\n",
]
for c in CAND:
    if c in s and ("88rpx" in ref or "babel" in ref):
        s = s.replace(c, "", 1)
        print(f"✓ 已移除（细节在 reference）: {c.strip()[:60]}")
    elif c in s:
        print(f"⚠ reference 里没找到对应细节，**保留**: {c.strip()[:60]}")

M.write_text(s, encoding="utf-8")
sz = M.stat().st_size
print(f"MEMORY.md -> {sz:,}/10240  余量 {10240-sz:+,}  {'OK' if sz <= 10240 else 'XX 超'}")
ms = M.read_text(encoding="utf-8")
refs = sorted(set(re.findall(r"`([\w\-./]+\.(?:md|py|npz|json))`", ms)))
dead = [p for p in refs if "/" in p and not ((WS / p).exists() or (ROOT / p).exists())]
print(f"引用 {len(refs)} 带路径死链 {len(dead)} {dead}")
for k in ("LRN-062", "罚因子伪影", "59,590"):
    print(f"  [{'OK' if k in ms else 'XX'}] MEMORY 含 {k}")
