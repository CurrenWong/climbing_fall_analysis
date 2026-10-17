"""记录一个关键发现：DRIVER=dispx 是 no-op（横向并未放开）。先记录再回答。"""
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")

# ── ANKLE_LIG_TODO.md：在 §3.1 表后补一段 ──
T = ROOT / "ANKLE_LIG_TODO.md"
s = T.read_text(encoding="utf-8")
anchor = "（`LOADDIR=z` 结果同族，无改善。）"
NOTE = anchor + """

### 3.1b 🔴 **更正（2026-10-04）：`DRIVER=dispx` 是 no-op，「横向自由」从未真正生效**

**证据**：M1(`all`+`dispx`)=**140 N** 与 M3(`all`+`disp`)=**140 N** 结果**完全相同**；
M2(`transverse`+`dispx`)=**60 N** 与 M4(`transverse`+`disp`)=**60 N** 也完全相同。

**根因**（`ankle_lig_feb.py` 880-895 行）：无论 `disp` 还是 `dispx`，
胫腓骨顶面都先写了一条 `zero displacement` 的 `_lateral` BC，
把**非轴向 DOF 全部钉死为 0**；`dispx` 只是**少写**轴向那一条。
⇒ 两种模式下顶面**实际上都被全约束**，`dispx` 的"横向自由"注释是**错的**。

**⇒ 后果**：**"放开横向 / 允许外翻+内旋"这个实验从未真正跑过。**
M1~M4 实际只扫了 `FIX_TALUS`（140 vs 60 N），不是顶面约束。

**真正能放开顶面的只有 `DRIVER=force` + `TOPFIX=none`**（888 行 `continue`），
但它**因刚体模态失败**（见 `FE_PIPELINE_TODO.md:41`）——大概率是
**初始接触未闭合时胫骨竖向不约束 ⇒ 奇异**。

**⇒ 正确的下一步**（`FE_PIPELINE_TODO.md` 实验 1 的可执行版）：
1. 修 `dispx`：`_lateral` BC 应跳过，只保留**最小**约束消刚体模态
   （如仅约束**单个**节点横向，或加弱弹簧）
2. 或走 `force` 路径：先规定一小段位移闭合间隙，再切 600 N 面力
3. 然后看四条韧带排序是否翻成 **吴恺序（胫跟＞胫舟＞跟腓＞下胫腓）**"""
assert anchor in s, "anchor not found"
s = s.replace(anchor, NOTE, 1)
T.write_text(s, encoding="utf-8")
print(f"ANKLE_LIG_TODO.md -> {T.stat().st_size:,} bytes")

# ── LEARNINGS.md：追加 LRN ──
L = ROOT / ".learnings/LEARNINGS.md"
l = L.read_text(encoding="utf-8")
if "LRN-20261004-038" not in l:
    l = l.rstrip() + """

## LRN-20261004-038: 「开关跑过了」≠「开关生效了」——先验证开关是否 no-op

- **类别**: `correction` / `knowledge_gap`
- **场景**: 踝关节 FE，`ankle_lig_feb.py` 有 `DRIVER=disp|dispx` 两个模式，
  代码注释写着 `dispx` 会让"横向自由：交给韧带 + 接触约束"。
  我据此认为"允许外翻/内旋的实验跑过了"。
- **实际**: `dispx` 是 **no-op**。两种模式都先写 `_lateral` BC 把非轴向 DOF 钉为 0，
  `dispx` 只少写轴向那一条 ⇒ 顶面实际都被全约束。
  **判据**：M1==M3（140 N）、M2==M4（60 N）**逐位相同**。
- **教训**: 当一个参数"理论上该有效果"但实测"两组结果逐位相同"时，
  **第一反应应该是"开关没生效"，而不是"这个因素不影响结果"。**
  逐位相同（数字+排序全同）是 no-op 的强指纹；真有物理效果通常会有微小差异。
- **行动**: ① 以后新增开关，必须用「开/关两组结果**必须有差异**」做冒烟测试
  ② 排查前先读 BC 代码，不信任注释 ③ 已回写 `ANKLE_LIG_TODO.md` §3.1b
- **关联**: `[LRN-20261003-027]`（位移控制/力控制病态）、`FE_PIPELINE_TODO.md:41`
"""
    L.write_text(l, encoding="utf-8")
    print(f"LEARNINGS.md -> {L.stat().st_size:,} bytes（+LRN-038）")
else:
    print("LEARNINGS.md 已有 LRN-038")
