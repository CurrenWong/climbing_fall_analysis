"""收口③：更新 FE_PIPELINE 迁移指针 + 全量校验。"""
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"
PIPE = os.path.join(R, "FE_PIPELINE_TODO.md")

t = io.open(PIPE, encoding="utf-8").read()
OLD = re.search(r"> ⚠️ \*\*本节的「踝子线」部分已迁移（2026-10-06）\*\*.*?"
                r"（加骨 von Mises 读数 → 骨强度判据 → 扫姿势出裕度表）。\n", t, re.S)
NEW = """> ⚠️ **踝子线已迁移（2026-10-06 收口）**：完整、自包含的状态与下一步见
> **`ANKLE_FIX_TODO.md` 顶部的「🧭 接手入口」**（12 行读完可接手）。三条必知：
> 1. **`dispx` 下 `DISP_MM>0` = 压实**（`[ERR-20261006-003]`）；零回归 `+1.6mm ⇒ 592.6020 N`。
> 2. 🎯 **相对结论可用、绝对数值不可用**：姿势排序经**两套软骨网格交叉验证**一致
>    ⇒「内翻 ⇒ 腓骨/外踝风险最高」可用；**绝对刚度带 2~4× 不确定度**，来源=三条建模伪影
>    （**缺距下关节 44.8% / 无摩擦切向滑移 39.5% / tie 罚滑移 13.9%**）。
> 3. **`[LRN-080]`「韧带 77%」已撤回**（~70% 功在罚弹簧）。
> **下一步：① 接入距下关节（软骨已生成⇒纯接线）② 力控找根接 OpenSim。**
"""
if OLD:
    t = t[:OLD.start()] + NEW + t[OLD.end():]
    io.open(PIPE, "w", encoding="utf-8").write(t)
    print(f"✓ FE_PIPELINE 迁移指针已更新 -> "
          f"{len(io.open(PIPE, encoding='utf-8').read().encode('utf-8'))} B")
else:
    print("⚠ 未匹配到旧指针（可能已更新）")

print("\n" + "=" * 74)
print("压缩前收口校验")
print("=" * 74)
CHK = [(os.path.join(W, "MEMORY.md"), 10240),
       (os.path.join(W, "references", "operational-gotchas.md"), None),
       (os.path.join(W, "references", "ankle-fe-criteria.md"), None),
       (os.path.join(W, "memory", "2026-10-06.md"), None),
       (os.path.join(R, "ANKLE_FIX_TODO.md"), None),
       (os.path.join(R, "FE_PIPELINE_TODO.md"), None),
       (os.path.join(R, "docs", "总体进度.md"), None),
       (os.path.join(R, "docs", "踝FE可用模型使用说明.md"), None),
       (os.path.join(R, ".learnings", "LEARNINGS.md"), None),
       (os.path.join(R, ".learnings", "ERRORS.md"), None)]
ok = True
for p, cap in CHK:
    n = len(io.open(p, encoding="utf-8").read().encode("utf-8"))
    f = ""
    if cap:
        f = " ✅" if n <= cap else f" ❌ 超 {n-cap}"
        ok &= n <= cap
    print(f"  {os.path.basename(p):<32} {n:>8,} B{f}")
for f in (".learnings/LEARNINGS.md", ".learnings/ERRORS.md"):
    x = io.open(os.path.join(R, f), encoding="utf-8").read()
    ids = re.findall(r"^#+\s*\[?([A-Z]{3}-\d{8}-\d{3})", x, re.M)
    dup = sorted({i for i in ids if ids.count(i) > 1})
    print(f"  {os.path.basename(f)}: {len(ids)} 条，标题级重号 {dup or '无 ✓'}")
    ok &= not dup
print("\n  --- 关键内容 ---")
need = [(os.path.join(R, "ANKLE_FIX_TODO.md"),
         ["🧭 接手入口（2026-10-06 收口", "相对可用 / 绝对不可用",
          "缺距下关节接触 44.8%", "韧带承担 77%」已撤回", "接入距下关节"]),
        (os.path.join(W, "MEMORY.md"),
         ["相对结论可用、绝对数值不可用", "三条建模伪影", "撤回 `[LRN-080]`",
          "缺距下关节接触 44.8%", "cart_m1_fix"]),
        (os.path.join(R, "FE_PIPELINE_TODO.md"),
         ["踝子线已迁移（2026-10-06 收口）", "两套软骨网格交叉验证"])]
for p, ks in need:
    tt = io.open(p, encoding="utf-8").read()
    miss = [k for k in ks if k not in tt]
    print(f"  {os.path.basename(p):<32} 缺失 {miss or '无 ✓'}")
    ok &= not miss
# 过时内容复查
stale = [(os.path.join(W, "MEMORY.md"), "韧带 77% / 接触 28%（最新）"),
         (os.path.join(W, "MEMORY.md"), "刚度 ~2× 未归因")]
for p, s in stale:
    tt = io.open(p, encoding="utf-8").read()
    print(f"  过时句 {s[:30]:<34} {'⚠ 仍在' if s in tt else '已清 ✓'}")
    ok &= s not in tt
print("\n结论：", "PASS ✅" if ok else "有问题 ⚠")
