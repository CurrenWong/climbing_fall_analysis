"""最终只读校验：新文档 + 链接 + MEMORY。"""
import re
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
NEW = ROOT / "docs/总体进度.md"

print("① 新文档")
print(f"   [{'OK' if NEW.exists() else 'MISS'}] {NEW.stat().st_size:,} B  docs/总体进度.md")
t = NEW.read_text(encoding="utf-8")
for k in ("线 1", "线 2", "子线", "P3", "V5", "A6", "LRN-051", "总体进度"):
    print(f"     [{'OK' if k in t else 'XX'}] {k}")

print("\n② git 可追踪性")
for f in ("docs/总体进度.md", "docs/方案更新_v2.md", "FE_PIPELINE_TODO.md"):
    r = subprocess.run(["git", "check-ignore", str(ROOT / f)], cwd=ROOT,
                       capture_output=True, text=True)
    tag = "被忽略 → 不进 git" if r.returncode == 0 else "可追踪 → 进 git"
    print(f"   {f:32s} {tag}")
r = subprocess.run(["git", "status", "--porcelain", "docs/总体进度.md"],
                   cwd=ROOT, capture_output=True, text=True)
print(f"   git status: {r.stdout.strip() or '(无变化)'}")

print("\n③ 链接")
for f in ("README.md", "FE_PIPELINE_TODO.md"):
    s = (ROOT / f).read_text(encoding="utf-8")
    print(f"   [{'OK' if '总体进度' in s else 'XX'}] {f}")

print("\n④ MEMORY.md")
M = WS / "MEMORY.md"
sz = M.stat().st_size
s = M.read_text(encoding="utf-8")
print(f"   disk = {sz:,} / 10240  余量 {10240 - sz:+,}  {'OK' if sz <= 10240 else '超限'}")
print(f"   [{'OK' if 'docs/总体进度.md' in s else 'XX'}] 总览指针")
refs = sorted(set(re.findall(r"`([\w\-./]+\.(?:md|py|npz|json|feb))`", s)))
dead = [p for p in refs
        if not ((WS / p).exists() or (ROOT / p).exists())
        and not p.startswith(("memory/", "references/", "digest/", "docs/"))
        and not list(ROOT.glob(f"**/{Path(p).name}"))]
print(f"   引用 {len(refs)} / 死链 {len(dead)} {dead}")

print("\n⑤ 本轮全部产物")
for k, p in [("docs/总体进度.md", NEW),
             ("README.md", ROOT / "README.md"),
             ("FE_PIPELINE_TODO.md", ROOT / "FE_PIPELINE_TODO.md"),
             ("ANKLE_FIX_TODO.md", ROOT / "ANKLE_FIX_TODO.md"),
             ("CONV_TODO.md", ROOT / "CONV_TODO.md"),
             ("cart_refine.py", ROOT / "scripts/ankle_fe/cart_refine.py"),
             ("LEARNINGS.md", ROOT / ".learnings/LEARNINGS.md"),
             ("MEMORY.md", M)]:
    print(f"   [{'OK ' if p.exists() else 'MISS'}] {p.stat().st_size if p.exists() else 0:>9,}  {k}")
lt = (ROOT / ".learnings/LEARNINGS.md").read_text(encoding="utf-8")
print(f"   LRN-051 {'OK' if 'LRN-20261004-051' in lt else 'XX'}；"
      f"LRN-050 更正声明 {'OK' if '本条已更正' in lt else 'XX'}")
