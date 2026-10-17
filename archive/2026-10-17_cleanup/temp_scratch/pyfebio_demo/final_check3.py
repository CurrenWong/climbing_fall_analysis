"""压缩上下文前的最终全量校验（纯 ASCII 输出，避免 cmd 吞行）。"""
import re
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")

print("=" * 78)
print("A. 文档清单与大小")
print("=" * 78)
DOCS = {
    "FE_PIPELINE_TODO.md (接手入口)": ROOT / "FE_PIPELINE_TODO.md",
    "ANKLE_LIG_TODO.md (韧带全部)": ROOT / "ANKLE_LIG_TODO.md",
    "docs/踝关节FE建模交接.md": ROOT / "docs/踝关节FE建模交接.md",
    ".learnings/LEARNINGS.md": ROOT / ".learnings/LEARNINGS.md",
    ".learnings/ERRORS.md": ROOT / ".learnings/ERRORS.md",
    "MEMORY.md": WS / "MEMORY.md",
    "memory/2026-10-03/ankle-fe-modeling.md":
        WS / "memory/2026-10-03/ankle-fe-modeling.md",
    "memory/2026-10-04/ankle-fe-solve-and-handover.md":
        WS / "memory/2026-10-04/ankle-fe-solve-and-handover.md",
    "temp/thums/ligament_results.txt": ROOT / "temp/thums/ligament_results.txt",
    "scripts/ankle_fe/ligaments.py": ROOT / "scripts/ankle_fe/ligaments.py",
    "scripts/ankle_fe/ankle_lig_feb.py": ROOT / "scripts/ankle_fe/ankle_lig_feb.py",
    "temp/pyfebio_demo/read_lig.py": ROOT / "temp/pyfebio_demo/read_lig.py",
}
missing = []
for k, p in DOCS.items():
    ok = p.exists()
    if not ok:
        missing.append(k)
    print(f"  [{'OK ' if ok else 'MISS'}] {p.stat().st_size if ok else 0:>9,}  {k}")

print()
print("=" * 78)
print("B. 关键内容锚点")
print("=" * 78)
mem = (WS / "MEMORY.md").read_text(encoding="utf-8")
fp = (ROOT / "FE_PIPELINE_TODO.md").read_text(encoding="utf-8")
al = (ROOT / "ANKLE_LIG_TODO.md").read_text(encoding="utf-8")
lrn = (ROOT / ".learnings/LEARNINGS.md").read_text(encoding="utf-8")
err = (ROOT / ".learnings/ERRORS.md").read_text(encoding="utf-8")
memnote = (WS / "memory/2026-10-04/ankle-fe-solve-and-handover.md"
           ).read_text(encoding="utf-8")

CHECKS = [
    ("MEMORY.md <= 10KB", len(mem.encode("utf-8")) <= 10240),
    ("MEMORY 有韧带轮条目", "韧带轮" in mem and "ANKLE_LIG_TODO" in mem),
    ("MEMORY 记负结果", "排序复现=负结果" in mem or "排序=负结果" in mem),
    ("接手入口 存在", "接手入口" in fp),
    ("入口 已更新韧带状态", "TRUNC_MM=60" in fp or "TRUNC_MM" in fp),
    ("入口 记 5.4 s", "5.4" in fp),
    ("ANKLE_LIG 有 4 韧带表", "81101203" in al and "81101101" in al),
    ("ANKLE_LIG 有 4 条件表", "Model initialization failed" in al),
    ("ANKLE_LIG 有负结果节", "负结果" in al),
    ("ANKLE_LIG 有下一步（图4）", "图 4" in al or "图4" in al),
    ("LRN-021 (补舟骨) 在档", "LRN-20261004-021" in lrn),
    ("LRN-027 (截断提速) 在档", "LRN-20261004-027" in lrn),
    ("LRN-028 (读迭代日志) 在档", "LRN-20261004-028" in lrn),
    ("LRN-029 (指标定义) 在档", "LRN-20261004-029" in lrn),
    ("ERR-004 (求解器配置) 在档", "ERR-20261004-004" in err),
    ("ERR-005 (pyfebio 补丁) 在档", "ERR-20261004-005" in err),
    ("当日记忆含韧带节", "加 8 组踝韧带" in memnote),
    ("当日记忆含 TRUNC_MM", "TRUNC_MM=60" in memnote),
]
for k, v in CHECKS:
    print(f"  [{'OK' if v else 'XX'}] {k}")

print()
print("=" * 78)
print("C. MEMORY.md 本地路径死链检查")
print("=" * 78)
paths = sorted(set(re.findall(r"`([\w\-./]+\.(?:md|py|npz|json|feb))`", mem)))
dead = [p for p in paths
        if not (WS / p).exists() and not (ROOT / p).exists()
        and not p.startswith("memory/2026-10-03/")
        and not p.startswith("memory/2026-10-04/")]
for p in paths:
    ok = ((WS / p).exists() or (ROOT / p).exists()
          or p.startswith(("memory/", "references/", "digest/")))
    print(f"  [{'OK ' if ok else 'DEAD'}] {p}")
print(f"  -> {len(paths)} 个引用，{len(dead)} 个可疑")

print()
allok = (not missing and all(v for _, v in CHECKS))
print(f"ALL PASS = {allok}")
print(f"MEMORY.md = {len(mem.encode('utf-8')):,} / 10240 bytes")
