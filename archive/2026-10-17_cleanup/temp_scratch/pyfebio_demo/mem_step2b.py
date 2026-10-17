"""MEMORY 再削 70+ B：gotcha #20 细节收进其指向文件（先验文件存在）。"""
from pathlib import Path
import re

WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
ROOT = Path(r"D:\Project\climbing_fall_analysis")
M = WS / "MEMORY.md"

# 先验证被指向的文件真的存在（搬不是删）
p20 = "digest/procedure/git-history-secret-purge.md"
found = [(b / p20) for b in (WS, ROOT) if (b / p20).exists()]
print(f"被指向文件存在: {[str(x) for x in found] or 'NO — 不能搬！'}")
if not found:
    raise SystemExit("目标文件不存在，拒绝删细节")

s = M.read_text(encoding="utf-8")
OLD = ("| 20 | **Git 历史敏感值清除** → `digest/procedure/git-history-secret-purge.md`："
       "`git bundle` 备份 → `git-filter-repo --replace-text` → fetch + `--force` |")
NEW = "| 20 | **Git 历史敏感值清除** → `digest/procedure/git-history-secret-purge.md` |"
if OLD in s:
    s = s.replace(OLD, NEW, 1)
    print("✓ gotcha #20 已瘦身")
else:
    print("⚠ 未匹配 #20")

M.write_text(s, encoding="utf-8")
sz = M.stat().st_size
print(f"MEMORY.md -> {sz:,}/10240  余量 {10240-sz:+,}  {'OK' if sz <= 10240 else 'XX 超'}")
# 死链复检
ms = M.read_text(encoding="utf-8")
refs = sorted(set(re.findall(r"`([\w\-./]+\.(?:md|py|npz|json))`", ms)))
dead = [p for p in refs if "/" in p and not ((WS / p).exists() or (ROOT / p).exists())]
print(f"引用 {len(refs)} 带路径死链 {len(dead)} {dead}")
