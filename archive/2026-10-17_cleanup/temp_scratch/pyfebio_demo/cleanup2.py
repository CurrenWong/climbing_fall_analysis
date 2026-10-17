"""精确删除 3 文件的双横幅 + 修 ANKLE_FIX_TODO 正文的最后一处相反论断。"""
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"
for rel in ("ANKLE_GOAL_TODO.md", "ANKLE_600N_TODO.md", "TIE_PORT_TODO.md"):
    p = os.path.join(R, rel)
    ls = io.open(p, encoding="utf-8").read().splitlines()
    i0 = next((i for i, l in enumerate(ls) if "2026-10-05 重大更正" in l), None)
    if i0 is None:
        print(f"· {rel} 无旧横幅")
        continue
    # 旧横幅到它自己的 '> ---' 及其后的空引用行为止
    j = i0
    while j < len(ls):
        if ls[j].strip() == "> ---":
            j += 1
            while j < len(ls) and ls[j].strip() in (">", ""):
                j += 1
            break
        j += 1
    removed = j - i0
    ls = ls[:i0] + ls[j:]
    io.open(p, "w", encoding="utf-8").write("\n".join(ls) + "\n")
    print(f"✓ {rel} 删除旧横幅 {removed} 行 -> {len(io.open(p, encoding='utf-8').read().encode('utf-8'))} B")

# ANKLE_FIX_TODO 正文的历史断定
P = os.path.join(R, "ANKLE_FIX_TODO.md")
t = io.open(P, encoding="utf-8").read()
OLD = "**⇒ 历史记录的 187.77 / 187.90 / 202.29 / 209.5589 / 592.61 N 全是拉伸值。**"
NEW = ("**⇒ ~~历史记录的 187.77 / 187.90 / 202.29 / 209.5589 / 592.61 N 全是拉伸值。~~ "
       "⚠️ 2026-10-06 再更正：这句**错了**（`[ERR-20261006-003]`）—— 那些是**压实**值，"
       "本条下面的 −132.56 / −608.82 N 才是拉伸。以顶部注记为准。**")
if OLD in t:
    t = t.replace(OLD, NEW)
    io.open(P, "w", encoding="utf-8").write(t)
    print(f"✓ ANKLE_FIX_TODO 正文矛盾已标出 -> {len(t.encode('utf-8'))} B")
else:
    print("⚠ ANKLE_FIX_TODO 目标句未命中")

# ── 复查 ──
print("\n" + "=" * 80)
print("最终复查：全部 8 个文件的错误方向论断")
print("=" * 80)
PAT = re.compile(r"(DISP_MM<0.{0,25}压缩|全是拉伸值|全是\*\*拉伸|LOAD_SIGN=\+1\*\* \| \*\*必开|"
                 r"DISP_MM=-1\.0.{0,12}压缩侧)")
bad = 0
for rel in ("docs/踝FE可用模型使用说明.md", "ANKLE_FIX_TODO.md", "docs/总体进度.md",
            "FE_PIPELINE_TODO.md", "ANKLE_GOAL_TODO.md", "ANKLE_600N_TODO.md",
            "TIE_PORT_TODO.md", "ANKLE_LIG_TODO.md"):
    p = os.path.join(R, rel)
    if not os.path.exists(p):
        continue
    t2 = io.open(p, encoding="utf-8").read()
    hits = [(i + 1, l.strip()) for i, l in enumerate(t2.splitlines()) if PAT.search(l)]
    print(f"  {rel:<38} {'✓ 干净' if not hits else str(len(hits)) + ' 处'}")
    for i, l in hits[:4]:
        print(f"       {i:>4}: {l[:130]}")
    bad += len(hits)
print(f"\n剩余疑似：{bad} 处")
