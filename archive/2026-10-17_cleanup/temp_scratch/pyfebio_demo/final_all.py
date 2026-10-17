"""记录"中文字符串混 ASCII 引号"坑 + 全量最终校验。"""
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"
GOT = os.path.join(W, "references", "operational-gotchas.md")
LEARN = os.path.join(R, ".learnings", "LEARNINGS.md")

# ── 记坑 ──
ADD = """

## 🐍 写脚本时的高频自伤（本项目一轮内各犯 3 次）

1. **中文字符串里混入 ASCII `"`** ⇒ `SyntaxError`。
   本会话在 3 个脚本上各犯一次（都是在写"引号引用的短语"时）。
   **⇒ 中文语境里一律用 `「」` 或 `『』`，外层用单引号。**
2. **cmd 里用多行 `python -c`** ⇒ `IndentationError`（同样 3 次）。
   **⇒ 多行一律写临时 `.py` 文件。**
3. **检查器的判据字符串写错** ⇒ 假阴性/假阳性而不自知：
   本轮实例 —— 只取脚本 stdout 时 `N O R M A L`（带空格）不出现 ⇒ 把正常收敛报成 ERROR；
   横幅检查用 `"必开" in text` ⇒ 把 `ALIGN_AXIS` 的"必开"误报为符号错误。
   **⇒ 判据上线前必须用一个"已知答案"的样例验一遍（见 `[LRN-20261006-083]` 规则 2）。**
"""
g = io.open(GOT, encoding="utf-8").read()
if "写脚本时的高频自伤" in g:
    print("· gotchas 已有")
else:
    io.open(GOT, "a", encoding="utf-8").write(ADD)
    print(f"operational-gotchas.md +{len(ADD)} 字符")

# ── 全量校验 ──
print("\n" + "=" * 74)
print("压缩前最终校验")
print("=" * 74)
import ast
ok = True
CHK = [(os.path.join(W, "MEMORY.md"), 10240), (os.path.join(W, "SOUL.md"), None),
       (os.path.join(W, "references", "operational-gotchas.md"), None),
       (os.path.join(W, "references", "ankle-fe-criteria.md"), None),
       (os.path.join(W, "memory", "2026-10-06.md"), None),
       (os.path.join(R, "docs", "踝FE可用模型使用说明.md"), None),
       (os.path.join(R, "ANKLE_FIX_TODO.md"), None),
       (os.path.join(R, "FE_PIPELINE_TODO.md"), None),
       (os.path.join(R, "docs", "总体进度.md"), None),
       (os.path.join(R, ".learnings", "LEARNINGS.md"), None),
       (os.path.join(R, ".learnings", "ERRORS.md"), None)]
for p, cap in CHK:
    if not os.path.exists(p):
        print(f"  ✗ 缺失 {p}")
        ok = False
        continue
    n = len(io.open(p, encoding="utf-8").read().encode("utf-8"))
    f = ""
    if cap:
        f = " ✅" if n <= cap else f" ❌ 超 {n-cap}"
        ok &= n <= cap
    print(f"  {os.path.basename(p):<32} {n:>8} B{f}")
for s in ("ankle_lig_feb.py", "ligaments.py", "check_c2.py", "lig_strain.py",
          "acceptance.py"):
    ast.parse(io.open(os.path.join(R, "scripts", "ankle_fe", s),
                      encoding="utf-8").read())
print("  AST ✓ ×5")
for f in (".learnings/LEARNINGS.md", ".learnings/ERRORS.md"):
    x = io.open(os.path.join(R, f), encoding="utf-8").read()
    ids = re.findall(r"^#+\s*\[?([A-Z]{3}-\d{8}-\d{3})", x, re.M)
    dup = sorted({i for i in ids if ids.count(i) > 1})
    print(f"  {os.path.basename(f)}: {len(ids)} 条，重号 {dup or '无 ✓'}")
    ok &= not dup
# 关键内容
need = [(os.path.join(R, "ANKLE_FIX_TODO.md"), ["恢复被误作废", "2.96 MPa", "287.4 mm²",
                                                "骨折判据", "P1"]),
        (os.path.join(R, "FE_PIPELINE_TODO.md"), ["本节的「踝子线」部分已迁移", "双靶已 PASS"]),
        (os.path.join(R, "docs", "总体进度.md"), ["骨折判据", "已建（首版）", "软骨网格健康度"]),
        (os.path.join(W, "MEMORY.md"), ["DISP_MM>0", "符号真相", "592.6"]),
        (os.path.join(R, ".learnings", "LEARNINGS.md"), ["旧配置的数字"])]
for p, ks in need:
    t = io.open(p, encoding="utf-8").read()
    miss = [k for k in ks if k not in t]
    print(f"  {os.path.basename(p):<32} 缺失 {miss or '无 ✓'}")
    ok &= not miss
# 旧错句复查
for rel in ("ANKLE_FIX_TODO.md", "ANKLE_GOAL_TODO.md", "ANKLE_600N_TODO.md",
            "TIE_PORT_TODO.md", "ANKLE_LIG_TODO.md", "FE_PIPELINE_TODO.md",
            "docs/总体进度.md", "docs/踝FE可用模型使用说明.md"):
    t = io.open(os.path.join(R, rel), encoding="utf-8").read()
    n_old = t.count("加载方向符号一直是反的")
    print(f"  旧错句 {rel:<34} {n_old} {'✓' if n_old == 0 else '⚠'}")
    ok &= n_old == 0
print("\n结论：", "PASS ✅" if ok else "有问题 ⚠")
