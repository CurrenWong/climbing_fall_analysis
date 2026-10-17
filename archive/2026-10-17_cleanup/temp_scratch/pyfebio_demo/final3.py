"""补日记（C3）+ 最终校验。"""
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"
DAY = os.path.join(W, "memory", "2026-10-06.md")

ADD = """
## 🎯 C3 也落地了：定量剂量-响应 + 阈值

**② 外翻不收敛查清了**：**6 组杠杆矩阵全部无效**（步长×10 / 落位限 0.1 / **落位限 0** /
`DET_MIN=1e-2` / 载荷减半），**只取决于角度**且失败在 **t→0（载荷≈0）**
⇒ **初态几何不兼容**：近退化软骨单元（`tibiotalar_talus` penta6 `det(J) min=2.39e-3`，
≈1000× 偏小）被压翻。**⇒ 属重分网/几何问题，不是求解器参数问题。**

**可用范围实测 `IE ∈ [−3°, +5°]`**（−3.5/−4/−5、+6/+8 全 ✗）

**③ 判据双侧 PASS —— "外侧−内侧"严格单调反对称**：

| IE | −3 | −2 | −1 | +1 | +2 | +3 | +4 | +5 |
|---|---|---|---|---|---|---|---|---|
| 外−内 | **−6.76** | −4.83 | −2.88 | +1.20 | +3.03 | +5.03 | +7.05 | **+9.07** |

**C3 定义**：`risk = ε_total(c) / ε_failure`，`ε_failure` 取自已有模块
（`src/climbing/coupling/ankle_ligament.py`，有测试断言）：ATFL 14%（11–15%）、CFL 13%（9–13%）
- CFL 风险比：+3° → 0.63（保守 0.90）、+4° → 0.77（保守 **1.11**）、+5° → 0.91（保守 **1.31**）
- **⇒ CFL 越阈点 +3°~+4°（保守下限）**；ATFL 本范围不越阈（≤0.36）

**⚠️ 绝对越阈点未定稿**（两个未定因素）：① CFL 建模 34.33 vs 解剖 ~28.5（稀释 ~20%）、
ATFL 36.99 vs 15~20（稀释 ~2×）⇒ 真值更高、越阈点前移；② 姿势预应变是模型化的
（零应变态定义在中性姿势）。**⇒ 现只可当"相对排序/相对风险"用。**

### 本轮完整落盘
- LRN-20261006-074~078（LEARNINGS）、ERR-20261006-001/002（ERRORS）
- `references/ankle-fe-criteria.md`（C2 三锚定 / 三口径实测 / **C3 剂量-响应** / 限制 / 工具）
- 使用说明 §4/§5 更新；ANKLE_FIX_TODO 追加第二轮
- 新工具：`check_c2.py`、`lig_strain.py`、`NO_RUN=1`

### 下一步（按价值）
1. **修软骨楔形单元退化**（detJ 2.4e-3）⇒ 打开外翻侧 −3.5~−5°
2. **修韧带插入点几何**（ATFL 36.99 → 15~20 mm）⇒ C3 绝对越阈点定稿
3. 把 C3 接到 OpenSim 的姿势输入（用户原目标）
"""
d = io.open(DAY, encoding="utf-8").read()
if "C3 也落地了" in d:
    print("· 日记已有该节")
else:
    io.open(DAY, "a", encoding="utf-8").write(ADD)
    print(f"✓ 日记已追加 +{len(ADD)} 字符")

print("\n" + "=" * 72)
print("最终校验（口径：解码后字节数）")
print("=" * 72)
ok = True
CHK = [(os.path.join(W, "MEMORY.md"), 10240),
       (os.path.join(W, "memory", "2026-10-06.md"), None),
       (os.path.join(W, "references", "ankle-fe-criteria.md"), None),
       (os.path.join(R, "docs", "踝FE可用模型使用说明.md"), None),
       (os.path.join(R, "ANKLE_FIX_TODO.md"), None),
       (os.path.join(R, ".learnings", "LEARNINGS.md"), None),
       (os.path.join(R, ".learnings", "ERRORS.md"), None),
       (os.path.join(R, "scripts", "ankle_fe", "check_c2.py"), None),
       (os.path.join(R, "scripts", "ankle_fe", "lig_strain.py"), None)]
for p, cap in CHK:
    t = io.open(p, encoding="utf-8").read()
    n = len(t.encode("utf-8"))
    f = ""
    if cap:
        f = " ✅" if n <= cap else f" ❌ 超 {n-cap}"
        ok &= n <= cap
    nm = os.path.relpath(p, R) if p.startswith(R) else os.path.relpath(p, W)
    print(f"  {nm:<46} {n:>8} B{f}")

import ast
for s in ("check_c2.py", "lig_strain.py", "ankle_lig_feb.py", "acceptance.py"):
    ast.parse(io.open(os.path.join(R, "scripts", "ankle_fe", s),
                      encoding="utf-8").read())
print("  AST ✓ ×4")

for f in (".learnings/LEARNINGS.md", ".learnings/ERRORS.md"):
    t = io.open(os.path.join(R, f), encoding="utf-8").read()
    ids = re.findall(r"^#+\s*\[?([A-Z]{3}-\d{8}-\d{3})", t, re.M)
    dup = sorted({i for i in ids if ids.count(i) > 1})
    print(f"  {os.path.basename(f)}: {len(ids)} 条，重号 {dup or '无 ✓'}")
    ok &= not dup

d2 = io.open(os.path.join(R, "docs", "踝FE可用模型使用说明.md"),
             encoding="utf-8").read()
need = ["C2 重定义", "+11.80%", "lig_strain.py", "三个口径", "内外翻不对称",
        "36.99 mm", "0.000°", "目前只有内翻侧可用", "CFL 越阈点"]
miss = [k for k in need if k not in d2]
print(f"  使用说明关键内容缺失：{miss or '无 ✓'}")
ok &= not miss

rr = io.open(os.path.join(W, "references", "ankle-fe-criteria.md"),
             encoding="utf-8").read()
needr = ["C2a-1", "三个口径", "+11.79", "−6.76", "越阈点"]
missr = [k for k in needr if k not in rr]
print(f"  references 关键内容缺失：{missr or '无 ✓'}")
ok &= not missr

print("\n结论：", "PASS ✅" if ok else "有问题 ⚠")
