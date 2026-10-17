"""把本轮结果追加到 ANKLE_FIX_TODO.md 与 daily note。"""
import io
import os

R = r"D:\Project\climbing_fall_analysis"
ADD = """

## 2026-10-05 多步落位（`[LRN-070]`）—— step1 打通，瓶颈转向几何

**已实现** `MULTISTEP=1`（默认关；零回归已验证：TIE+CART_NEW 仍 **592.6020 N NORMAL**）。
`<Step>`：step1 位移落位 `SEAT_MM` → step2 力控制 `FORCE_N`。**step1 收敛成功（t=1）**。

**修 2 个真 bug**：
1. `symmetric_stiffness` **必须在 `<solver>` 内**（放 Control 直下报 `unrecognized tag`）。
   刚体连接件的硬要求，缺它 ⇒ step2 立刻 47 负 J；修后刚度矩阵 nnz 258222→507067。
2. **`<traction>` 不写 `lc` ⇒ FEBio 把满载荷在第一个增量全部加上**（实测 step2 首个增量 48 负 J）。
   修法：加第二条曲线 `(0,0)(1,0)(2,1)` + `lc="2"` ⇒ **负 J 归零**。

**🔴 剩余瓶颈 = 几何/网格，不是边界条件**：
- 软骨只能承受 ~0.8 mm 压缩就翻负（`dispx` 刚体距骨：0.5mm=351.51N NORMAL、0.9mm 失败）
- t=0 关节间隙中位 **~1.0 mm** ⇒ **压不闭** ⇒ 间隙未闭时胫骨轴向只靠韧带 ⇒ step2 ~12 N 即不收敛
- **刚体距骨把可用压缩范围砍半**（0.8 vs 可变形 1.6 mm），压缩全集中在软骨
- **下一步**：让关节在 t=0 闭合（软骨厚度填满关节间隙）—— Anderson 正常踝本身即接触态
"""
P = os.path.join(R, "ANKLE_FIX_TODO.md")
s = io.open(P, encoding="utf-8").read()
io.open(P, "a", encoding="utf-8").write(ADD)
print(f"ANKLE_FIX_TODO.md: {len(s)} -> {len(s)+len(ADD)}")

D = os.path.join(R, "memory", "2026-10-05.md")
DAILY = """
## 踝线多步落位（LRN-070）
- **MULTISTEP=1 落地**：step1 位移落位收敛成功（t=1）→ step2 力控制
- 修 2 bug：`symmetric_stiffness` 须在 `<solver>` 内；`<traction>` 无 `lc` 会用满载荷
- **瓶颈转为几何**：软骨 ~0.8 mm 就翻负，而初始间隙 ~1.0 mm ⇒ 压不闭
- 零回归：TIE+CART_NEW = 592.6020 N NORMAL
"""
old = io.open(D, encoding="utf-8").read() if os.path.exists(D) else "# 2026-10-05\n"
io.open(D, "a", encoding="utf-8").write(DAILY)
print(f"daily: {len(old)} -> {len(old)+len(DAILY)}")
for f in ("ANKLE_FIX_TODO.md", "memory/2026-10-05.md",
          ".learnings/LEARNINGS.md", "scripts/ankle_fe/ankle_lig_feb.py"):
    fp = os.path.join(R, f)
    txt = io.open(fp, encoding="utf-8").read()
    print(f"  {f:<34} {os.path.getsize(fp):>8} B  LRN-070×{txt.count('LRN-20261005-070')}"
          f"  MULTISTEP×{txt.count('MULTISTEP')}")
