"""MEMORY.md 追加 tied-elastic 结果（注意 10 KB 硬上限）。"""
from pathlib import Path

M = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md")
LIMIT = 10240
s = M.read_text(encoding="utf-8")

OLD = ("- ⭐ **根因=接触区欠分辨**：仅 **73 节点**(~2.5mm)；加密到 692 后软骨体积涨\n"
       "  **2.3×** ⇒ 粗网格低估面积")
NEW = ("- ⭐ **根因=接触区欠分辨**：仅 **73 节点**(~2.5mm)；加密到 692 后软骨体积涨\n"
       "  **2.3×** ⇒ 粗网格低估面积\n"
       "- ✅ **tied-elastic 验证通过**（连非共形网格，`cart_refine.py` + `TIE=1`）：\n"
       "  L=0 罚法**单调收敛**到共形解 436.817 N → penalty 1/2/**4** = 404.7/420.9/**432.5**\n"
       "  （**−0.99%**）；8+ 首步崩。坑：tie 面必须是**真实 facet**（tet4 只能三角形）；\n"
       "  **secondary 有洞 ⇒ 伪装成「硬 2.2×」**（支反力随 penalty 上升背离真值）；\n"
       "  ⚠️ 软骨面内加密仍被**母单元扭曲**挡住（hex8 加密恒 30% 子单元翻负）")
if "tied-elastic 验证通过" not in s:
    assert OLD in s, "锚点未找到"
    s = s.replace(OLD, NEW, 1)
    M.write_text(s, encoding="utf-8")
sz = len(s.encode())
print(f"disk = {sz:,} / {LIMIT}   余量 {LIMIT - sz:+,}")
if sz > LIMIT:
    s = s.replace("（连非共形网格，`cart_refine.py` + `TIE=1`）", "（`cart_refine.py`+`TIE=1`）")
    s = s.replace("  L=0 罚法**单调收敛**到共形解 436.817 N → penalty 1/2/**4** = 404.7/420.9/**432.5**\n  （**−0.99%**）；8+ 首步崩。坑：tie 面必须是**真实 facet**（tet4 只能三角形）；\n  **secondary 有洞 ⇒ 伪装成「硬 2.2×」**（支反力随 penalty 上升背离真值）；\n  ⚠️ 软骨面内加密仍被**母单元扭曲**挡住（hex8 加密恒 30% 子单元翻负）",
                  "  L=0 罚法**单调收敛**到共形解 436.817 N → penalty 1/2/**4**=404.7/420.9/**432.5**\n"
                  "  （**−0.99%**）。坑：tie 面须是**真实 facet**（tet4 只三角形）；**secondary 有洞\n"
                  "  ⇒ 伪装成「硬 2.2×」**；⚠️ 面内加密仍被**母单元扭曲**挡住（hex8 恒 30% 翻负）")
    M.write_text(s, encoding="utf-8")
    sz = len(s.encode())
    print(f"裁剪 -> {sz:,}   余量 {LIMIT - sz:+,}")
print(f"最终 {sz:,} / {LIMIT}  {'OK' if sz <= LIMIT else '超限'}")
