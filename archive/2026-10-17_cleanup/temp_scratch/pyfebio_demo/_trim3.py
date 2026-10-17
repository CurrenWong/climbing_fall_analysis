"""最后一轮瘦身：last_updated 压缩 + 索引描述精简 + SecLog 表合并。"""
from pathlib import Path

P = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md")
t = P.read_text(encoding="utf-8")
before = len(t.encode("utf-8"))
done = []


def sub(old, new, tag):
    global t
    assert t.count(old) == 1, f"[{tag}] 命中 {t.count(old)}"
    t = t.replace(old, new)
    done.append(tag)


# ① last_updated 压到两行
i = t.index("**最后更新**:")
j = t.index("**整理者**:")
t = t[:i] + (
    "**最后更新**: 2026-10-04　"
    "**① MEMORY.md 整理** 19.5 KB → 9.7 KB：SecLog/CloudBase 12 条 gotcha 详情搬进新 "
    "`references/seclog-miniprogram-gotchas.md`（主表改「症状→一句话根因」速查）；"
    "gotcha 20/21 收成指针；**修掉一处原有死链**（markitdown-ocr 配置实际在 skill 里）。　"
    "**② 踝关节 FE 第二轮**：「提高软骨网格质量」= **负面结果**"
    "（`|V|min` 提升 42 倍但收敛更差 ⇒ 瓶颈是接触，`LRN-020`）；"
    "subtalar 关节**网格建好但接触面非 2-流形**（`LRN-019`）；"
    "单关节回归 rc=0/16.6 s/**436.8 N**，软骨 p99 2.58–4.27 MPa。\n"
) + t[j:]
done.append("last_updated")

# ② references 索引：描述精简
sub("| `references/seclog-miniprogram-gotchas.md` | **SecLog 小程序 / CloudBase 全部 gotcha 详版**（2026-10-04 从主 MEMORY 迁出）|",
    "| `references/seclog-miniprogram-gotchas.md` | **SecLog 小程序 / CloudBase 全部 gotcha 详版** |", "idx1")
sub("| `references/cloudbase-ai.md` | **CloudBase AI 多模态调用**：cloud.ai() 正确路径、SDK 版本坑、ApiKey 频控误判 |",
    "| `references/cloudbase-ai.md` | CloudBase AI 多模态：cloud.ai() 正确路径、SDK 版本坑、频控误判 |", "idx2")
sub("| `references/tools-cheatsheet.md` | **65 个 skill 速查**（33 agent 内置 / 6 active 脚本 / 26 全局 npx）|",
    "| `references/tools-cheatsheet.md` | **65 个 skill 速查**（33 内置 / 6 脚本 / 26 npx）|", "idx3")
sub("| `references/miniprogram-ci.md` | 微信小程序编译 + 预览二维码（私钥路径 + ES2020 `?.` 兼容坑 + agent-ui 补丁模式）|",
    "| `references/miniprogram-ci.md` | 小程序编译 + 预览二维码（私钥、`?.` 兼容、agent-ui 补丁）|", "idx4")
sub("| `references/fe-modeling-gotchas.md` | **FE 建模陷阱**：接触四铁律 / tet 翻转奇偶性 / FEBio `.feb` 语法 / pyfebio 坑 / 二分方法论 |",
    "| `references/fe-modeling-gotchas.md` | **FE 陷阱**：接触四铁律 / tet 奇偶性 / `.feb` 语法 / pyfebio 坑 / 二分法 |", "idx5")
sub("| `references/operational-gotchas.md` | 进程/超时/删除/Self-Improvement/定时任务 |",
    "| `references/operational-gotchas.md` | 进程/超时/删除/Self-Improvement/定时 |", "idx6")

# ③ SecLog 表：行级合并（避开全角引号）
lines = t.split("\n")
i_list = next(i for i, l in enumerate(lines) if "month:'this'" in l and l.startswith("|"))
lines[i_list] = ("| 日期跨月错位（月初统计=0 / 本月混入上月） | 库存 UTC + 云函数跑 UTC ⇒ "
                 "月边界按 `TZ_OFFSET=8h` 构造；`list` 传 `month:'this'` 而非 `days` 滚动窗口 |")
i_mon = next(i for i, l in enumerate(lines) if "TZ_OFFSET=8h` 构造" in l
             and "跨月错位" not in l)
lines[i_mon] = ""
done.append("tbl1")
t = "\n".join(lines)

# ④ 头部
sub("> **硬上限 10 KB**。超了就是没守住纪律（2026-10-04 曾在 19.5 KB，整理回 9.8 KB）。",
    "> **硬上限 10 KB**（2026-10-04 曾到 19.5 KB，已整理）。", "head")

P.write_text(t, encoding="utf-8")
after = len(t.encode("utf-8"))
print(f"已应用 {len(done)} 处: {done}")
print(f"MEMORY.md {before:,} -> {after:,} bytes "
      f"({'OK, 余量 ' + str(10240 - after) if after <= 10240 else 'FAIL'})")
