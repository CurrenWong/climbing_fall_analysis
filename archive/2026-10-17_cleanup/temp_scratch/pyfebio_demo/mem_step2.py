"""MEMORY.md：写入本轮关键结论（dispx+transverse 稳定配置 / Loads bug），并等量削减。"""
from pathlib import Path

M = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md")
s = M.read_text(encoding="utf-8")

OLD = """- 🎯 **子线下一步（已定）**：目标=**接 OpenSim + 判不同姿势**。卡口=**力控制病态**
  （距骨全固定⇒刚体模态；`FIX_TALUS=all`⇒双端夹持⇒踝不能内外翻）。解=**Anderson §2.2
  弱弹簧 500/50 N/mm + 解除内外翻内旋**；验收=**姿势排序对罚因子鲁棒**。
  ⚠️ 韧带线结案是**循环论证**（ΔL≈0 是双端夹持强制的）⇒ 待重开。`[LRN-057]`"""
NEW = """- 🎯 **子线下一步（已定）**：目标=**接 OpenSim + 判不同姿势**。⚠️「力控制病态」是**误诊**：
  真因=① **`<Loads>` 恒为空**（`continue` 后死代码 ⇒ `force` 载荷被静默丢，已修）② **关节初始张开**
  （0.6 N 也失败 463 s）。⇒ **稳定姿势配置=`DRIVER=dispx`+`FIX_TALUS=transverse`**（NORMAL/5s/0 负J，
  与 `all` 不同 90 vs 162 N）；**读数出口=`LOG_DATA`**；验收=**姿势排序对罚因子鲁棒**。`[LRN-057~059]`
  ⚠️ 韧带线结案是**循环论证**（ΔL≈0 由双端夹持强制）⇒ 待重开"""

TRIMS = [
    ("""### 文档 / 内容
- 公众号发布 / SVG→PNG 无白边 → `references/wechat-publish.md`
- PDF/DOCX 转换 → `references/docling.md`（10 页以上走批量脚本）
- OCR + 文本提取 → `skills/markitdown-ocr/`
- 小程序编译 / 预览二维码 → `references/miniprogram-ci.md`""",
     """### 文档 / 内容
- 公众号 / SVG→PNG → `references/wechat-publish.md`；PDF/DOCX → `references/docling.md`（10 页+ 走批量）
- OCR → `skills/markitdown-ocr/`；小程序编译/预览码 → `references/miniprogram-ci.md`"""),
    ("""### 其他
- MCP 死代码 → 自接管：模板 `references/mcp-bypass.md`（Tavily / Flint 已确认同款）

- 定时任务 → `copaw cron`（skill）；Server 酱推送 → `references/operational-gotchas.md`""",
     """### 其他
- MCP 死代码 → 自接管：模板 `references/mcp-bypass.md`（Tavily / Flint 同款）
- 定时任务 → `copaw cron`（skill）；Server 酱推送 → `references/operational-gotchas.md`"""),
]
if OLD in s:
    s = s.replace(OLD, NEW, 1)
    print("✓ 踝段已更新")
else:
    print("⚠ 未匹配踝段")
for a, b in TRIMS:
    if a in s:
        s = s.replace(a, b, 1)
    else:
        print(f"⚠ 未匹配削减: {a[:34]}...")
M.write_text(s, encoding="utf-8")
sz = M.stat().st_size
t = M.read_text(encoding="utf-8")
print(f"MEMORY.md -> {sz:,}/10240  余量 {10240-sz:+,}  {'OK' if sz <= 10240 else 'XX 超'}")
for k in ("dispx", "Loads", "LRN-057~059", "循环论证"):
    print(f"  [{'OK' if k in t else 'XX'}] {k}")
