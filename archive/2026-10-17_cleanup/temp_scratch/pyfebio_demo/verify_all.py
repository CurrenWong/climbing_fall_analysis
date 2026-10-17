"""补回紧凑的「详见」指针 + 最终全量校验。"""
import re
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
M = WS / "MEMORY.md"
LIMIT = 10240

s = M.read_text(encoding="utf-8")
# 在踝关节段末尾（下一节标题前）插入紧凑指针
anchor = "|精度 **10~15%**"
if "详见 `memory/2026-10-03" not in s:
    i = s.find(anchor)
    j = s.find("\n", i)
    s = s[:j] + "\n- **详见**：`memory/2026-10-03/ankle-fe-modeling.md`§六–§十" + s[j:]
    M.write_text(s, encoding="utf-8")
sz = M.stat().st_size
print(f"补指针后 disk = {sz:,} / {LIMIT}   余量 {LIMIT - sz:+,}")
if sz > LIMIT:
    s = s.replace("§六–§十", "", 1)
    M.write_text(s, encoding="utf-8")
    sz = M.stat().st_size
    print(f"  再裁 -> {sz:,}   余量 {LIMIT - sz:+,}")

print("\n" + "=" * 74)
print("① 文档清单")
DOCS = {
    "docs/Anderson双论文_关键信息.md": ROOT / "docs/Anderson双论文_关键信息.md",
    "docs/吴恺论文_重解析关键信息.md": ROOT / "docs/吴恺论文_重解析关键信息.md",
    "docs/踝关节FE建模交接.md": ROOT / "docs/踝关节FE建模交接.md",
    "ANKLE_FIX_TODO.md": ROOT / "ANKLE_FIX_TODO.md",
    "CONV_TODO.md": ROOT / "CONV_TODO.md",
    "ANKLE_600N_TODO.md": ROOT / "ANKLE_600N_TODO.md",
    "ANKLE_LIG_TODO.md": ROOT / "ANKLE_LIG_TODO.md",
    "FE_PIPELINE_TODO.md": ROOT / "FE_PIPELINE_TODO.md",
    "FE_PIPELINE_TODO.md 入口指针": None,
    ".learnings/LEARNINGS.md": ROOT / ".learnings/LEARNINGS.md",
    ".learnings/ERRORS.md": ROOT / ".learnings/ERRORS.md",
    "scripts/ankle_fe/ankle_lig_feb.py": ROOT / "scripts/ankle_fe/ankle_lig_feb.py",
    "scripts/ankle_fe/build_feb.py": ROOT / "scripts/ankle_fe/build_feb.py",
    "scripts/ankle_fe/cartilage.py": ROOT / "scripts/ankle_fe/cartilage.py",
    "temp/thums/thums_lowerlimb_r1.npz": ROOT / "temp/thums/thums_lowerlimb_r1.npz",
    "temp/thums/r1/cartilage.json": ROOT / "temp/thums/r1/cartilage.json",
}
miss = []
for k, p in DOCS.items():
    if p is None:
        continue
    ok = p.exists()
    if not ok:
        miss.append(k)
    print(f"  [{'OK ' if ok else 'MISS'}] {p.stat().st_size if ok else 0:>9,}  {k}")

print("\n② 本轮 LRN 全在档")
lrn = (ROOT / ".learnings/LEARNINGS.md").read_text(encoding="utf-8")
n_ok = 0
for n in (39, 40, 41, 42, 43, 44, 45):
    k = f"LRN-20261004-0{n}"
    ok = k in lrn
    n_ok += ok
    print(f"  [{'OK' if ok else 'XX'}] {k}")

print("\n③ MEMORY.md 死链")
paths = sorted(set(re.findall(r"`([\w\-./]+\.(?:md|py|npz|json|feb|k))`", s)))
dead = []
for p in paths:
    if (WS / p).exists() or (ROOT / p).exists() or \
       p.startswith(("memory/", "references/", "digest/", "docs/")):
        continue
    hit = list(ROOT.glob(f"**/{Path(p).name}"))
    if hit:
        continue
    dead.append(p)
print(f"  {len(paths)} 引用 / {len(dead)} 死链  {dead}")

print("\n④ 关键结论是否在 MEMORY 里")
for nm, k in [("刚度 2 倍", "259.3"), ("根因=接触区欠分辨", "73 节点"),
              ("Anderson 真值", "all-hex"), ("韧带线降级", "Tochigi"),
              ("CONV 指针", "CONV_TODO"), ("修正计划", "ANKLE_FIX_TODO")]:
    print(f"  [{'OK' if k in s else 'XX'}] {nm}")

print(f"\nALL = {'PASS' if not miss and not dead and n_ok == 7 else 'CHECK'}")
print(f"MEMORY.md = {sz:,} / {LIMIT}")
