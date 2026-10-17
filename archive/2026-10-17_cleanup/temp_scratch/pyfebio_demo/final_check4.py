"""压缩上下文前的最终全量校验（纯 ASCII 输出）。"""
import re
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")

print("=" * 78)
print("A. 文档清单")
print("=" * 78)
DOCS = {
    "CONV_TODO.md (收敛研究/关键问题)": ROOT / "CONV_TODO.md",
    "ANKLE_600N_TODO.md (600N)": ROOT / "ANKLE_600N_TODO.md",
    "ANKLE_LIG_TODO.md (韧带)": ROOT / "ANKLE_LIG_TODO.md",
    "FE_PIPELINE_TODO.md (接手入口)": ROOT / "FE_PIPELINE_TODO.md",
    "docs/吴恺论文_重解析关键信息.md": ROOT / "docs/吴恺论文_重解析关键信息.md",
    ".learnings/LEARNINGS.md": ROOT / ".learnings/LEARNINGS.md",
    ".learnings/ERRORS.md": ROOT / ".learnings/ERRORS.md",
    "MEMORY.md": WS / "MEMORY.md",
    "memory/…/ankle-fe-solve-and-handover.md":
        WS / "memory/2026-10-04/ankle-fe-solve-and-handover.md",
    "fepatch.py (复用补丁)": ROOT / "temp/pyfebio_demo/fepatch.py",
}
missing = []
for k, p in DOCS.items():
    ok = p.exists()
    if not ok:
        missing.append(k)
    print(f"  [{'OK ' if ok else 'MISS'}] {(p.stat().st_size if ok else 0):>9,}  {k}")

print()
print("=" * 78)
print("B. 关键锚点")
print("=" * 78)
conv = (ROOT / "CONV_TODO.md").read_text(encoding="utf-8")
fp = (ROOT / "FE_PIPELINE_TODO.md").read_text(encoding="utf-8")
mem = (WS / "MEMORY.md").read_text(encoding="utf-8")
note = (WS / "memory/2026-10-04/ankle-fe-solve-and-handover.md"
        ).read_text(encoding="utf-8")
lrn = (ROOT / ".learnings/LEARNINGS.md").read_text(encoding="utf-8")
p600 = (ROOT / "ANKLE_600N_TODO.md").read_text(encoding="utf-8")

CHECKS = [
    ("CONV 有「关键问题」节", "关键问题" in conv),
    ("CONV 三层齐全", all(x in conv for x in ("第一层", "第二层", "第三层"))),
    ("CONV 有最尖锐候选(接触面定向)", "接触面定向" in conv),
    ("CONV 有方法论问题(拆开排查)", "拆开" in conv),
    ("CONV 记 259.3/116.9", "259.3" in conv and "116.9" in conv),
    ("CONV 记五假设", all(x in conv for x in ("①", "②", "③", "④", "⑤"))),
    ("FE_PIPELINE 入口有刚度警告", "2 倍" in fp and "CONV_TODO" in fp),
    ("MEMORY 有 2 倍警告", "2 倍" in mem and "CONV_TODO" in mem),
    ("MEMORY 记双靶可用", "双靶可用" in mem),
    ("MEMORY 记载荷-位移不可信", "载荷-位移不可信" in mem),
    ("600N TODO 记双靶结果", "2.96" in p600 and "287.4" in p600),
    ("记忆有 §八 收敛研究", "§八" in note or "八、收敛研究" in note),
    ("记忆记五假设排除", "五个假设" in note or "①" in note),
    ("LRN-033/034/035/036/037 在档",
     all(f"LRN-20261004-0{n}" in lrn for n in (33, 34, 35, 36, 37))),
]
for k, v in CHECKS:
    print(f"  [{'OK' if v else 'XX'}] {k}")

print()
print("=" * 78)
print("C. MEMORY.md 路径死链")
print("=" * 78)
paths = sorted(set(re.findall(r"`([\w\-./]+\.(?:md|py|npz|json|feb))`", mem)))
dead = 0
for p in paths:
    ok = ((WS / p).exists() or (ROOT / p).exists()
          or p.startswith(("memory/", "references/", "digest/")))
    if not ok:
        dead += 1
        print(f"  [DEAD] {p}")
print(f"  {len(paths)} 个引用，{dead} 个死链")

print()
allok = (not missing and all(v for _, v in CHECKS) and dead == 0)
print(f"ALL PASS = {allok}")
print(f"MEMORY.md = {len(mem.encode('utf-8')):,} / 10240 bytes")
