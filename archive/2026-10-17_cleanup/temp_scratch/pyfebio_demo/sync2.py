"""MEMORY.md 裁到磁盘≤10240；FE_PIPELINE_TODO.md 表格加 tied-elastic 行。"""
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
M = WS / "MEMORY.md"
L = 10240
d = lambda p: p.stat().st_size  # noqa: E731

# ══════ ① FE_PIPELINE_TODO.md：在「三件必须知道」表里加一行 ══════
F = ROOT / "FE_PIPELINE_TODO.md"
ft = F.read_text(encoding="utf-8")
ROW = ("| ✅ **新通路（2026-10-04）** | **`tied-elastic` 验证可用**（手册 §3.14.6，连接非共形网格）：\n"
       "L=0 罚法**单调从下方收敛**到共形解 436.817 N → penalty 1/2/**4** = 404.7/420.9/**432.5 N**\n"
       "（**−0.99%**），8+ 首步即崩。工具 `scripts/ankle_fe/cart_refine.py` + `build_feb.py` 的 `TIE=1`。\n"
       "⚠️ 但**软骨面内红加密本身仍被母单元扭曲挡住**（hex8 加密恒 30% 子单元翻负，`[LRN-050]`）\n"
       "⇒ **机制通、收益尚未兑现**。三个坑 → `ANKLE_FIX_TODO.md` P1d2a、`[LRN-046~049]` |")
if "tied-elastic" not in ft:
    anchor = "| 🚧 **已知阻塞** |"
    i = ft.find(anchor)
    j = ft.find("\n", i)
    ft = ft[:j + 1] + ROW + "\n" + ft[j + 1:]
    F.write_text(ft, encoding="utf-8")
print(f"FE_PIPELINE_TODO.md -> {d(F):,} B")

# ══════ ② MEMORY.md 裁到磁盘 ≤10240 ══════
s = M.read_text(encoding="utf-8")
print(f"MEMORY.md 起始 disk = {d(M):,} (余量 {L - d(M):+,})")
CUTS = [
    # 把上次误加长的「最后更新」前缀删掉
    ("**最后更新**: 2026-10-04　**⓪ tied-elastic 验证**（→`ANKLE_FIX_TODO.md` P1d2a）："
     "L=0 罚法单调收敛 404.7/420.9/**432.5 N**（共形 436.817，**−0.99%**），甜点 "
     "`penalty=4`；坑=LRN-046~050。　**① MEMORY.md 整理**",
     "**最后更新**: 2026-10-04　**⓪ tied-elastic 已验证**（→`ANKLE_FIX_TODO.md` P1d2a）"
     "：penalty=4 → 432.5 N（共形 436.8，−0.99%）。　**① MEMORY.md 整理**"),
    ("| 20 | **Git 历史敏感值清除** → `digest/procedure/git-history-secret-purge.md`："
     "`git bundle` 备份 → `git-filter-repo --replace-text` → fetch + `--force`。"
     "微信 App ID 不是密钥；filter-repo 会删 origin；`--force-with-lease` 必被拒 |",
     "| 20 | **Git 历史敏感值清除** → `digest/procedure/git-history-secret-purge.md`："
     "`git bundle` 备份 → `git-filter-repo --replace-text` → fetch + `--force` |"),
    ("| 21 | **FE 接触四铁律 + tet 翻转奇偶性 + `.feb` 语法** | `references/fe-modeling-gotchas.md` + `[LRN-20261003-016~020]`。**一句话**：",
     "| 21 | **FE 接触四铁律 + tet 翻转奇偶性 + `.feb` 语法** | `references/fe-modeling-gotchas.md`。**一句话**："),
    ("| 22 | **`build_feb.py` ≠ `ankle_lig_feb.py`**（前者无开关）；**cmd `set` 污染整条命令行** | `[LRN-042/045]` |",
     "| 22 | **`build_feb.py` ≠ `ankle_lig_feb.py`**（前者无开关）；**cmd `set` 污染整条命令行** | `[LRN-042/045]` |"),
]
for a, b in CUTS:
    if d(M) <= L:
        break
    if a in s and a != b:
        s = s.replace(a, b, 1)
        M.write_text(s, encoding="utf-8")
        print(f"  -> {d(M):,} (余量 {L - d(M):+,})")

MORE = [
    ("🔴 **根因=接触区欠分辨**", "⭐ **根因=接触区欠分辨**"),
    ("  软骨量级 · k=1000 N/mm 测 0.07~0.16 mm", "  软骨量级 · k=1000 测 0.07~0.16 mm"),
    ("- **`uv pip install gmsh` ≠ 官方 SDK**", "- **`uv pip install gmsh` ≠ 官方 SDK**"),
]
for a, b in MORE:
    if d(M) <= L:
        break
    if a in s and a != b:
        s = s.replace(a, b, 1)
        M.write_text(s, encoding="utf-8")
        print(f"  more -> {d(M):,}")

print(f"\nMEMORY.md disk = {d(M):,} / {L}  {'OK' if d(M) <= L else '仍超 ' + str(d(M) - L)}")
