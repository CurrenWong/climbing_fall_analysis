"""收官：压 维护规则 与 last_updated。"""
from pathlib import Path

P = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md")
t = P.read_text(encoding="utf-8")
before = len(t.encode("utf-8"))

# ① 维护规则整段重写（1,209 -> ~600）
i = t.index("## \U0001F504 维护规则")
j = t.index("**最后更新**:")
t = t[:i] + """## \U0001F504 维护规则

1. **MEMORY.md ≤ 10 KB 硬上限**
2. 新内容默认进 `references/`，主 MEMORY 只留指针；项目专属细节归对应 references
3. **搬，不是删** —— 迁走前先确认新文件里有全部信息
4. 整理完必验三项：**尺寸达标 + 无死链（每个路径都存在）+ 信息未丢**
5. 更新时同步改 `last_updated`；重大重构前先备份到 `memory_archive/`

""" + t[j:]

# ② last_updated 压到 ~380
i = t.index("**最后更新**:")
j = t.index("**整理者**:")
t = t[:i] + (
    "**最后更新**: 2026-10-04　"
    "**① MEMORY.md 整理** 19.5 KB → 9.6 KB：SecLog/CloudBase 12 条 gotcha 详情搬进新 "
    "`references/seclog-miniprogram-gotchas.md`（主表改「症状→一句话根因」速查）；"
    "gotcha 20/21 收成指针；修掉一处原有死链。　"
    "**② 踝关节 FE 第二轮**：「提高软骨网格质量」= **负面结果**"
    "（`|V|min` 提升 42 倍但收敛更差 ⇒ 瓶颈是接触，`LRN-020`）；"
    "subtalar 关节网格建好但**接触面非 2-流形**（`LRN-019`）；"
    "单关节回归 rc=0/16.6 s/**436.8 N**，软骨 p99 2.58–4.27 MPa。\n"
) + t[j:]

P.write_text(t, encoding="utf-8")
after = len(t.encode("utf-8"))
print(f"MEMORY.md {before:,} -> {after:,} bytes "
      f"({'OK, 余量 ' + str(10240 - after) if after <= 10240 else 'FAIL, 差 ' + str(after-10240)})")
