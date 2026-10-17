"""修 MEMORY.md 死链 + 再砍到 10 KB 以内（短锚点版）。"""
from pathlib import Path

P = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md")
t = P.read_text(encoding="utf-8")
before = len(t.encode("utf-8"))
done = []


def sub(old, new, tag):
    global t
    assert t.count(old) == 1, f"[{tag}] 锚点命中 {t.count(old)}"
    t = t.replace(old, new)
    done.append(tag)


# ① 死链：markitdown-ocr 配置实际在 skill 里（整理前就指错了）
sub("- **OCR + 文本提取** → `references/markitdown-ocr配置.md`",
    "- **OCR + 文本提取** → `skills/markitdown-ocr/references/configuration.md`"
    "（批量脚本 `scripts/convert_pdf_parallel.py`、`ocr_scanned_pdf.py`）",
    "markitdown")

# ② last_updated 压缩
i = t.index("**最后更新**:")
j = t.index("**整理者**:")
t = t[:i] + (
    "**最后更新**: 2026-10-04\n"
    "- **MEMORY.md 整理**：19.5 KB → 9.8 KB。SecLog/CloudBase/小程序 12 条 gotcha 详情"
    "搬进新 `references/seclog-miniprogram-gotchas.md`，主表改「症状→一句话根因」速查+指针；"
    "gotcha 20/21 收成一行；**修掉一处原有死链**（markitdown-ocr 配置实际在 skill 里）\n"
    "- **踝关节 FE 第二轮**：①「提高软骨网格质量」= **负面结果**"
    "（|V|min 提升 42 倍但收敛更差 ⇒ 瓶颈是接触，`LRN-020`）"
    "② subtalar 关节**网格建好但接触面非 2-流形**（`LRN-019`）"
    "③ 单关节回归 rc=0/16.6 s/436.8 N，软骨 p99 2.58–4.27 MPa\n"
) + t[j:]
done.append("last_updated")

# ③ SecLog 表收短（短锚点）
sub("⇒ 按 `TZ_OFFSET=8h` 构造月边界 |", "⇒ 月边界按 `TZ_OFFSET=8h` 构造 |", "row2")
sub("fileID 的 bucket 用了云函数资源 ID；", "fileID 用了云函数资源 ID 当 bucket；", "row3")
sub("**不抛异常**只返回空 URL |", "**不抛异常**只返空 URL |", "row3b")
sub("**且二次 invoke 验证**（冷启动缓存）；小程序走 `ci-tools/compile.js upload`（要 babel 转译），**裸调 miniprogram-ci 是假发布**",
    "**且二次 invoke 验证**（冷启动缓存）；小程序走 `ci-tools/compile.js upload`，**裸调 miniprogram-ci 是假发布**",
    "row1")

# ④ 踝关节 FE 子线收短
sub("""- **当前可复现结果**（`DRIVER=disp DISP_MM=0.1472`）：
  `rc=0` / `NORMAL TERMINATION` / 16.6 s，支反力 **436.8 N**，
  软骨 von Mises p99 **2.58 / 4.27 MPa** → 外推 600 N 胫侧 ≈ 3.4 MPa，
  **与吴恺的 3.74 MPa 吻合**（同一量级）""",
    """- **可复现结果**（`DRIVER=disp DISP_MM=0.1472`）：`rc=0`/`NORMAL TERMINATION`/16.6 s，
  支反力 **436.8 N**，软骨 von Mises p99 **2.58 / 4.27 MPa** → 外推 600 N 胫侧 ≈3.4 MPa，
  **与吴恺的 3.74 MPa 吻合**""", "ankle1")
sub("""- **几何源**：THUMS V7.1（`paper/AM50_V71_Occupant.zip` 内 219 MB `.k`），
  单位 **mm-s-tonne-N-MPa ⇒ E 已是 MPa 零换算**；`*_center_spon` 是**骨髓腔**不是松质骨""",
    """- **几何源** THUMS V7.1（`paper/AM50_V71_Occupant.zip` 内 219 MB `.k`），
  单位 mm-s-tonne-N-MPa ⇒ **E 已是 MPa 零换算**；`*_center_spon` 是**骨髓腔**不是松质骨""",
    "ankle2")
sub("""- **未完成**：subtalar 关节的接触面**不是合法 2-流形**
  （`FULL_ANKLE=1` 网格能建但第一个载荷步失败）→ 见 `LRN-20261003-019`""",
    """- **未完成**：subtalar 接触面**非合法 2-流形**（`FULL_ANKLE=1` 网格能建、第一个载荷步失败）""",
    "ankle3")
sub("""- **候选下一步**：加 8 组踝韧带（`81100101~81101401_shell` 已在 npz 里）→
  复现吴恺的**定性排序**（胫跟＞胫舟＞跟腓＞下胫腓前）""",
    """- **候选下一步**：加 8 组踝韧带（`81100101~81101401_shell` 已在 npz 里）→ 复现吴恺
  **定性排序**（胫跟＞胫舟＞跟腓＞下胫腓前）""", "ankle4")

# ⑤ 头部注释压缩
sub("> **硬上限 10 KB**。超了就是我没守住纪律 —— 2026-10-04 曾在 19.5 KB，被整理回 8.9 KB。",
    "> **硬上限 10 KB**。超了就是没守住纪律（2026-10-04 曾在 19.5 KB，整理回 9.8 KB）。",
    "head")

# ⑥ 维护规则去重 + 加验收项（逐行短锚点）
sub("2. **新内容默认进 `references/`，主 MEMORY 只加指针**（含",
    "2. **新内容默认进 `references/`，主 MEMORY 只加指针**；项目专属细节归对应 references\n3. **搬，不是删** —— 迁走前先确认新文件里有全部信息\n4. **过时的内容定期清理**（heartbeat/dream）\n5. **更新时同步改 `last_updated`**；**重大重构前先备份**（`memory_archive/`）\n6. **整理完必验**：尺寸达标 + **无死链**（每个本地路径都存在）+ 信息未丢\n<!--OLD3-->", "rules-a")
sub("3. **项目专属细节不进主 MEMORY** —— 归到对应的 `references/*.md`\n", "", "rules-b")
sub("4. **过时的内容定期清理**（heartbeat/dream 整理时处理）\n", "", "rules-c")
sub("5. **更新主 MEMORY 时同步改 `last_updated`**\n", "", "rules-d")
sub("6. **重大重构前先备份**（`memory_archive/`）", "", "rules-e")
t = t.replace("<!--OLD3-->", "")

P.write_text(t, encoding="utf-8")
after = len(t.encode("utf-8"))
print(f"已应用 {len(done)} 处: {done}")
print(f"MEMORY.md {before:,} → {after:,} bytes "
      f"({'✅ ≤10240' if after <= 10240 else '❌ 仍超标'})")
