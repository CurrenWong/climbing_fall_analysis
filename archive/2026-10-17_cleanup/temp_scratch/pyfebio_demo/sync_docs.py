"""同步 tied-elastic 结果到文档；并把 MEMORY.md 按**磁盘口径**裁到 ≤10240。

★ 口径纠正：MEMORY.md 是 CRLF，`len(text.encode())` 比 `stat().st_size` 少
  每行 1 字节（145 行 ⇒ 少 145）。**必须用 stat().st_size 判定**，否则会误以为达标。
"""
import io
import re
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
L = 10240


def disk(p: Path) -> int:
    return p.stat().st_size


# ══════════ ① MEMORY.md 裁到磁盘 ≤10240 ══════════
M = WS / "MEMORY.md"
s = M.read_text(encoding="utf-8")
print(f"MEMORY.md 起始 disk = {disk(M):,} / {L}")

CUTS = [
    ("- ✅ **tied-elastic 验证通过**（`cart_refine.py`+`TIE=1`）:L=0 罚法**单调收敛**到共形解\n"
     "  436.817 N → penalty 1/2/**4**=404.7/420.9/**432.5**（**−0.99%**）。\n"
     "  坑：tie 面须是**真实 facet**（tet4 只三角）；**secondary 有洞 ⇒ 伪装成「硬 2.2×」**；\n"
     "  ⚠️ 面内加密仍被**母单元扭曲**挡住（hex8 恒 30% 翻负，`[LRN-046~050]`）",
     "- ✅ **tied-elastic 验证通过**（`cart_refine.py`+`TIE=1`）：L=0 罚法**单调收敛**到共形解\n"
     "  436.817 N → penalty 1/2/**4**=404.7/420.9/**432.5**（**−0.99%**）。坑：tie 面须是\n"
     "  **真实 facet**（tet4 只三角）；**secondary 有洞 ⇒ 伪装成「硬 2.2×」**；\n"
     "  ⚠️ 面内加密仍被**母单元扭曲**挡住（hex8 恒 30% 翻负，`[LRN-046~050]`）"),
    ("| 20 | **Git 历史敏感值清除**：`git bundle` 备份 → `git-filter-repo --replace-text` → fetch + `--force` 强推 | `digest/procedure/git-history-secret-purge.md`。要点：微信 App ID 不是密钥；filter-repo 会删 origin 要重加；`--force-with-lease` 必被拒；`git grep $(git rev-list --all)` 核验不可信 |",
     "| 20 | **Git 历史敏感值清除** → `digest/procedure/git-history-secret-purge.md`：`git bundle` 备份 → `git-filter-repo --replace-text` → fetch + `--force`。微信 App ID 不是密钥；filter-repo 会删 origin；`--force-with-lease` 必被拒 |"),
    ("| 19~20 |", "| 19~20 |"),
    ("**最后更新**: 2026-10-04　**① MEMORY.md 整理**",
     "**最后更新**: 2026-10-04　**⓪ tied-elastic 验证**（→`ANKLE_FIX_TODO.md` P1d2a）："
     "L=0 罚法单调收敛 404.7/420.9/**432.5 N**（共形 436.817，**−0.99%**），甜点 "
     "`penalty=4`；坑=LRN-046~050。　**① MEMORY.md 整理**"),
]
for a, b in CUTS:
    if a in s and a != b:
        s = s.replace(a, b, 1)
        M.write_text(s, encoding="utf-8")
        print(f"  -> disk {disk(M):,} (余量 {L - disk(M):+,})")

# 兜底：还超就逐条压
FALLBACK = [
    ("  软骨 E=12 ν=0.42、厚 1.5~1.7mm｜精度 **10~15%**",
     "  软骨 E=12 ν=0.42、厚 1.5~1.7mm｜精度 10~15%"),
    ("⇒ 我 ΔL≈0 **与实验一致**；吴恺 1.92~5.24mm 隐含 8~21% **不符实验**",
     "⇒ 我 ΔL≈0 **与实验一致**；吴恺 1.92~5.24mm 隐含 8~21% **不符实验**"),
    ("**同版本号可能是两份二进制 → 判据是实测 API 可调用性**",
     "**同版本号可能是两份二进制；判据=实测 API 可调用性**"),
]
for a, b in FALLBACK:
    if disk(M) <= L:
        break
    if a in s:
        s = s.replace(a, b, 1)
        M.write_text(s, encoding="utf-8")
        print(f"  fallback -> disk {disk(M):,}")

print(f"最终 MEMORY.md disk = {disk(M):,} / {L}  {'OK' if disk(M) <= L else '超限'}")

# ══════════ ② FE_PIPELINE_TODO.md：入口加 tied-elastic 状态 ══════════
F = ROOT / "FE_PIPELINE_TODO.md"
t = F.read_text(encoding="utf-8")
if "tied-elastic" not in t:
    # 在「🧭 接手入口」节的第一条状态行后插入
    m = re.search(r"(##[^\n]*🧭[^\n]*\n)", t)
    ins = ""
    if m:
        ins = ("\n> **最新（2026-10-04）**：`tied-elastic`（手册 §3.14.6）**已验证可用**——\n"
               "> L=0 罚法单调从下方收敛到共形解 436.817 N：penalty 1/2/**4** =\n"
               "> 404.7/420.9/**432.5 N（**−0.99%**），8+ 首步崩。\n"
               "> 工具：`scripts/ankle_fe/cart_refine.py` + `build_feb.py` 的 `TIE=1`；\n"
               "> 详情与三个坑 → `ANKLE_FIX_TODO.md` P1d2a，`[LRN-046~050]`。\n"
               "> ⚠️ **软骨面内加密本身仍被母单元扭曲挡住**（P1d2b）⇒ 目前**尚不能**\n"
               "> 用它在提升接触分辨率上兑现收益。\n")
        t = t[:m.end()] + ins + t[m.end():]
    else:
        t = ins + t
    F.write_text(t, encoding="utf-8")
print(f"\nFE_PIPELINE_TODO.md -> {disk(F):,}  (插入入口指针)")

# ══════════ ③ CONV_TODO.md：交叉引用 ══════════
C = ROOT / "CONV_TODO.md"
c = C.read_text(encoding="utf-8")
if "tied-elastic" not in c:
    c = c.rstrip() + (
        "\n\n---\n\n## 相关（2026-10-04）\n\n"
        "`tied-elastic` 已用同一套 600 N 模型验证可用（L=0 罚法单调收敛到共形解\n"
        "436.817 N，penalty=4 → 432.497 N，−0.99%）—— 见 `ANKLE_FIX_TODO.md` P1d2a\n"
        "与 `[LRN-046~050]`。它给「接触区欠分辨」提供了一条**不动骨网格**的加密通路，\n"
        "但**软骨面内加密本身**目前被母单元扭曲挡住（P1d2b），故尚未能兑现。\n")
    C.write_text(c, encoding="utf-8")
print(f"CONV_TODO.md -> {disk(C):,}  (加交叉引用)")
