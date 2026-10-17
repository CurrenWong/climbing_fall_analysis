"""追加 daily note（memory 目录在 agent workspace 下）。"""
import io
import os

D = r"C:\Users\curre\.copaw\workspaces\default\memory\2026-10-05.md"
DAILY = """
## 踝线多步落位（LRN-070）—— step1 打通，瓶颈转向几何
- **`MULTISTEP=1` 落地**（默认关）：step1 位移落位 → step2 力控制；**step1 收敛成功（t=1）**
- 修 2 个真 bug：
  1. `symmetric_stiffness` **必须在 `<solver>` 内**（放 Control 直下报 unrecognized tag）；
     刚体连接件硬要求，缺它 ⇒ step2 立刻 47 负 J；修后 nnz 258222→507067
  2. **`<traction>` 不写 `lc` ⇒ FEBio 把满载荷在第一个增量全加上**（step2 首个增量 48 负 J）；
     修法：第二条曲线 `(0,0)(1,0)(2,1)` + `lc="2"` ⇒ 负 J 归零
- **🔴 瓶颈转为几何**：软骨 ~0.8 mm 压缩就翻负，而 t=0 关节间隙中位 ~1.0 mm ⇒ **压不闭**
  ⇒ 间隙未闭时胫骨轴向只靠韧带 ⇒ step2 ~12 N 即不收敛；**刚体距骨把可用压缩砍半**（0.8 vs 1.6 mm）
- 零回归：TIE+CART_NEW = **592.6020 N NORMAL**
"""
old = io.open(D, encoding="utf-8").read() if os.path.exists(D) else "# 2026-10-05\n"
io.open(D, "a", encoding="utf-8").write(DAILY)
print(f"daily {len(old)} -> {len(old)+len(DAILY)}  ({D})")
