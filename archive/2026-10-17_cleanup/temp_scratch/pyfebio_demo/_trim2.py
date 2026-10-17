"""二轮瘦身 v3：全部用极短锚点。"""
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


sub("（批量脚本 `scripts/convert_pdf_parallel.py`、`ocr_scanned_pdf.py`）",
    "（批量脚本在 `skills/markitdown-ocr/scripts/`）", "ocr")

# gotcha 21：两处收短
sub("（我用", "（我用", "noop") if False else None
sub("改善了 42 倍却让收敛更差）|", "|", "g21b")
sub("**接触面必须与单元网格的真实面重合且是合法 2-流形**（否则接触完全不检测，"
    "但求解报 `rc=0` 无警告、**只是应力全 0**）",
    "**接触面要与网格真实面重合且是合法 2-流形**"
    "（否则接触不检测，却报 `rc=0` 无警告、**只是应力全 0**）", "g21a")

# 工具入口：逐条
sub("- **微信公众号发布** → `references/wechat-publish.md`\n"
    "- **PDF/DOCX 转换** → `references/docling.md`（10 页以上 → 批量脚本）\n"
    "- **OCR + 文本提取** → `skills/markitdown-ocr/references/configuration.md`"
    "（批量脚本在 `skills/markitdown-ocr/scripts/`）\n"
    "- **SVG → PNG 无白边** → `references/wechat-publish.md` §SVG\n"
    "- **微信小程序编译/预览二维码** → `references/miniprogram-ci.md`",
    "- 公众号发布 / SVG→PNG 无白边 → `references/wechat-publish.md`\n"
    "- PDF/DOCX 转换 → `references/docling.md`（10 页以上走批量脚本）\n"
    "- OCR + 文本提取 → `skills/markitdown-ocr/`\n"
    "- 小程序编译 / 预览二维码 → `references/miniprogram-ci.md`", "tools1")
sub("### MCP 死代码 → 自接管\n"
    "- 通用模板 → `references/mcp-bypass.md`；Tavily / Flint 是已确认两个同款问题",
    "### 其他\n"
    "- MCP 死代码 → 自接管：模板 `references/mcp-bypass.md`（Tavily / Flint 已确认同款）",
    "tools2")
sub("### 心跳 / 提醒 / 调度\n"
    "- 定时任务管理 → `copaw cron`（skill）；Server 酱推送 → "
    "`references/operational-gotchas.md` §必做",
    "- 定时任务 → `copaw cron`（skill）；Server 酱推送 → "
    "`references/operational-gotchas.md`", "tools3")

# 通用工具坑
sub("（无 `addDiscreteSurface`/`reclassify`/`geo.addTriangle`）⇒ STL→实体走不通。\n"
    "  要用官方 Windows64 SDK（`D:\\Program\\gmsh-sdk\\`，DLL 在 `lib/` 不在 `bin/`）。\n"
    "  **通用教训：同版本号可能是两份不同二进制 → 判断依据是实测 API 可调用性**",
    "（occ 绑定残缺 ⇒ STL→实体走不通）；官方 SDK 在 `D:\\Program\\gmsh-sdk\\`，\n"
    "  DLL 在 `lib/` 不在 `bin/`。**通用教训：同版本号可能是两份不同二进制 → "
    "判据是实测 API 可调用性**", "gmsh")
sub("**cmd.exe 下别用 `python -c \"多行\"`**（会 `IndentationError`）——\n"
    "  多行脚本一律写临时 `.py` 再跑，用完删\n"
    "- **cmd 下没有 `tail`**：用 `more` 或重定向到文件再 `read_file`\n"
    "- **edit 工具 payload 过大/含特殊符号时会失败，且可能\"部分应用\"留下损坏文件**\n"
    "  ⇒ 改完一律 `ast.parse()` 验语法，别只看工具返回值",
    "**cmd 下别用 `python -c \"多行\"`**（`IndentationError`）⇒ 写临时 `.py` 再跑；\n"
    "  **cmd 没有 `tail`** ⇒ 用 `more` 或重定向后 `read_file`\n"
    "- **edit 工具 payload 过大/含特殊符号会失败，且可能部分应用留下损坏文件**\n"
    "  ⇒ 改完一律 `ast.parse()` 验语法，别只看工具返回值", "cmd")

# 写作偏好：先删两条独立行，再合并中文优先行
sub("- **先行动后报告**：能跑就跑，跑出问题再回看\n", "", "st1")
sub("- **宁缺毋滥**：搜索类任务硬性 3-5 条，凑数扣分\n", "", "st2")
sub("- **中文优先**：标题/正文用中文，专有名词保留英文\n",
    "- **中文优先**（专有名词保留英文）；**先行动后报告**（能跑就跑，跑出问题再回看）；\n"
    "  **宁缺毋滥**（搜索硬性 3-5 条，凑数扣分）\n", "st3")

P.write_text(t, encoding="utf-8")
after = len(t.encode("utf-8"))
print(f"已应用 {len(done)} 处: {done}")
print(f"MEMORY.md {before:,} -> {after:,} bytes "
      f"({'OK, 余量 ' + str(10240 - after) if after <= 10240 else 'FAIL'})")
