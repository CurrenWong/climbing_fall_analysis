from pathlib import Path

M = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md")
L = 10240
s = M.read_text(encoding="utf-8")
CUTS = [
    ("- **🗺️总体进度**：`docs/总体进度.md`（两条线+子线一页总览）",
     "- **🗺️总览**：`docs/总体进度.md`（两条线+子线）"),
    ("- **入口**：`FE_PIPELINE_TODO.md`「🧭接手」｜`ANKLE_FIX_TODO.md`｜`docs/Anderson双论文_关键信息.md`",
     "- **入口**：`FE_PIPELINE_TODO.md`「🧭接手」｜`ANKLE_FIX_TODO.md`｜`docs/Anderson双论文_关键信息.md`"),
    ("**最后更新**: 2026-10-04　**⓪ tied-elastic 已验证**（→`ANKLE_FIX_TODO.md` P1d2a）"
     "：penalty=4 → 432.5 N（共形 436.8，−0.99%）。　**① MEMORY.md 整理**",
     "**最后更新**: 2026-10-04　**⓪ tied-elastic 已验证**（→`ANKLE_FIX_TODO.md` P1d2a）"
     "：penalty=4 → 432.5 N（−0.99%）。　**① MEMORY.md 整理**"),
    ("- **详见**：`memory/2026-10-03/ankle-fe-modeling.md`",
     "- **详见**：`memory/2026-10-03/ankle-fe-modeling.md`"),
    ("  ⚠️ 面内加密仍被**母单元扭曲**挡住（hex8 恒 30% 翻负，`[LRN-046~050]`）",
     "  ⚠️ 面内加密仍被**母单元扭曲**挡住（hex8 恒 28% 翻负，`[LRN-046~051]`）"),
]
for a, b in CUTS:
    if M.stat().st_size <= L:
        break
    if a in s and a != b:
        s = s.replace(a, b, 1)
        M.write_text(s, encoding="utf-8")
print(f"disk = {M.stat().st_size:,} / {L}  余量 {L - M.stat().st_size:+,}")
