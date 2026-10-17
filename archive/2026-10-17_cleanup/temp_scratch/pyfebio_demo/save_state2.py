"""追加本轮第二轮结论到 ANKLE_FIX_TODO + daily。"""
import io
import os

R = r"D:\Project\climbing_fall_analysis"
ADD = """

### 2026-10-05 第二轮（`[ERR-20261005-008]` + `[LRN-20261005-071]`）

**🐛 找到并修掉最重要的一个 bug：`lc` 挂错元素**
- 写成 `<traction lc="2">` ⇒ **FEBio 静默忽略** ⇒ **满载荷在第一增量瞬加**
- 正解（手册 §3.13.2.2 / p121）：`<scale lc="2">1.0</scale>` + `<traction>…</traction>`
- **⇒ 此前历史上所有 `DRIVER=force` 实验都没真正测过"斜坡力控制"**
- 诊断信号（值得记）：**失败点对总载荷量级不敏感**（100/600/5000 N 都在第一增量失败）
  ⇒ 必是"满值瞬加"，不是"斜坡"

**✅ 已验证正确（推翻上一轮的两条推测）**
- **多步机制正确**：`FORCE_N≈0` 时 NORMAL 跑完 43 个时间点，且 step2 顶面反力≈0
  ⇒ step1 的位移约束**确实**在 step2 被正确解除
- **关节 t=0 几乎已接触**（推翻"间隙 1 mm 压不闭"）：SEAT 0.02 mm 就有 31 N 反力；
  0.05→43.7、0.10→65.2、0.20→109.0、0.30→155.0、0.50→245.1 N（胫骨顶）
  ⇒ **不需要改几何**
- **诊断口径**：`reac_*.txt` 是每节点一行，须 `reshape(nstep,nnode,3).sum(1)`；第一列是节点 ID

**🟡 仍未通过**：`SEAT_HOLD=1`（step2 保持约束）把 step2 从 1.01818(10.9N) 推进到
**1.055(33.0N)** ⇒ "撤约束回弹"是真因素之一但非全部；~33 N 处仍 `Negative jacobian`。

**下一步（按价值排序）**
1. **按 Anderson 原设计把约束加在「距骨」而非胫骨顶** —— 距骨垂直位移当常驻运动学把手，
   胫骨顶从不被位移约束、只受力 ⇒ 无回弹，且忠于 Anderson §2.2
2. 打开 `ls_check_jacobians` 查 step1 收敛态是否本身已含翻负单元
3. 软骨网格加密（`cart_patch` m≥2 建网格即失败，待攻）
"""
P = os.path.join(R, "ANKLE_FIX_TODO.md")
s = io.open(P, encoding="utf-8").read()
io.open(P, "a", encoding="utf-8").write(ADD)
print(f"ANKLE_FIX_TODO.md {len(s)} -> {len(s)+len(ADD)}")

D = r"C:\Users\curre\.copaw\workspaces\default\memory\2026-10-05.md"
DAILY = """
### 踝线多步落位 · 第二轮
- **🐛 关键 bug**：FEBio 的 `<traction lc="...">` **被静默忽略** ⇒ 满载荷瞬加。
  正解是 `<scale lc="2">1.0</scale>`（手册 §3.13.2.2）。**历史上所有 force 实验都没测过真斜坡。**
  诊断信号：失败点对载荷量级不敏感（100/600/5000N 都首增量失败）
- ✅ 已验证：多步机制正确（FORCE_N≈0 时 NORMAL 跑完）；**关节 t=0 几乎已接触**
  （SEAT 0.02mm→31N、0.50mm→245N）⇒ **不需改几何**（推翻上一轮推测）
- 🟡 SEAT_HOLD=1 把 step2 从 10.9N 推到 33N ⇒ 回弹是真因素之一但非全部
- 下一步：**把约束加在距骨而非胫骨**（Anderson §2.2 原设计）
"""
o = io.open(D, encoding="utf-8").read()
io.open(D, "a", encoding="utf-8").write(DAILY)
print(f"daily {len(o)} -> {len(o)+len(DAILY)}")
