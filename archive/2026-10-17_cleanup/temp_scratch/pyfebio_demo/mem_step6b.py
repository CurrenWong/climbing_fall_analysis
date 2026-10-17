"""MEMORY 削 ~330 B：#21 细节搬进 reference（先验证 detail 在 reference 里）。"""
from pathlib import Path
import re

WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
ROOT = Path(r"D:\Project\climbing_fall_analysis")
M = WS / "MEMORY.md"
REF = WS / "references/fe-modeling-gotchas.md"

ref = REF.read_text(encoding="utf-8") if REF.exists() else ""
print(f"reference: {REF.exists()}  {len(ref):,} chars")
for key in ("刚体模态", "2-流形", "奇置换", "瓶颈"):
    print(f"  [{'OK' if key in ref else 'XX'}] reference 含 {key}")

s = M.read_text(encoding="utf-8")
OLD = ("| 21 | **FE 接触四铁律 + tet 翻转奇偶性 + `.feb` 语法** | `references/fe-modeling-gotchas.md`。"
       "**一句**：① 只靠接触约束的部件必有**刚体模态**⇒奇异⇒单元成批翻负 ⇒ 加载端**整个截面全约束** "
       "② **接触面要与网格真实面重合且是合法 2-流形**（否则静默不检测，只报 `rc=0`、应力全 0） "
       "③ **tet4 翻转必须奇置换**（整体倒序是偶置换 ⇒ det 符号不变！） "
       "④ **优化前先确认它是不是瓶颈**（曾误判\"网格质量\"，提升 42 倍反而更差）|")
NEW = ("| 21 | **FE 四铁律**（接触/翻转/`.feb` 语法）+ **\"某驱动不行先验产物再看报错\"** | "
       "`references/fe-modeling-gotchas.md`。① \"只靠接触⇒刚体模态\"是**不完整诊断**（`[LRN-058]`："
       "真因常是载荷被丢/初始未咬合）② 接触面须与真实网格面重合且合法 2-流形 ③ tet4 翻转须**奇置换** "
       "④ 优化前先确认它是瓶颈 |")
if OLD in s:
    s = s.replace(OLD, NEW, 1)
    print("✓ #21 已瘦身")
elif "| 21 |" in s:
    i = s.find("| 21 |")
    j = s.find("\n", i)
    s = s[:i] + NEW + s[j:]
    print("✓ #21 已替换（按行）")
else:
    print("⚠ 未找到 #21")
M.write_text(s, encoding="utf-8")
sz = M.stat().st_size
print(f"MEMORY.md -> {sz:,}/10240  余量 {10240-sz:+,}  {'OK' if sz <= 10240 else 'XX 超'}")
ms = M.read_text(encoding="utf-8")
refs = sorted(set(re.findall(r"`([\w\-./]+\.(?:md|py|npz|json))`", ms)))
dead = [p for p in refs if "/" in p and not ((WS / p).exists() or (ROOT / p).exists())]
print(f"引用 {len(refs)} 带路径死链 {len(dead)} {dead}")
for k in ("踝FE可用模型使用说明.md", "POSE_TILT", "LRN-063", "59,590"):
    print(f"  [{'OK' if k in ms else 'XX'}] MEMORY 含 {k}")
