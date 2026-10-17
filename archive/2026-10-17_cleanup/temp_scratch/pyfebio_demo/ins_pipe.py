"""在 FE_PIPELINE_TODO 的「接手入口」标题后插入迁移提示。"""
import io
import os

R = r"D:\Project\climbing_fall_analysis"
P = os.path.join(R, "FE_PIPELINE_TODO.md")
ls = io.open(P, encoding="utf-8").read().splitlines()
NOTE = [
    "",
    '> ⚠️ **本节的「踝子线」部分已迁移（2026-10-06）**：完整、自包含的状态与计划见',
    '> **`ANKLE_FIX_TODO.md` §第六轮**。三个必须知道的变化：',
    '> 1. **压实 = `DISP_MM > 0`**（`[ERR-20261006-003]`，三证）。**历史值 592.6/187.77/202.29 N 是压实、有效**；',
    '>    `DISP<0` 的 −132…−1017 N 才是拉伸（已作废）。600 N 在 `DISP_MM≈+1.9 mm`。',
    '> 2. **双靶已 PASS**：2.96 MPa(−21%) / 287.4 mm²(−14%) vs 吴恺 3.74 MPa / 335.5 mm²。',
    '> 3. **C2 PASS、C3 可绝对化**（力臂与解剖吻合：CFL 28.5/ATFL 21/PTFL 8 mm）。',
    '>',
    '> **最大缺口 = 骨折判据**（用户目标是「判不同姿势**骨折**风险」，而此前做的是韧带/扭伤）',
    '> ⇒ 计划 P1（加骨 von Mises 读数 → 骨强度判据 → 扫姿势出裕度表）。',
    "",
]
found = False
for i, l in enumerate(ls):
    if l.startswith("# 🧭 接手入口"):
        ls[i + 1:i + 1] = NOTE
        found = True
        break
if found:
    io.open(P, "w", encoding="utf-8").write("\n".join(ls) + "\n")
    print(f"✓ FE_PIPELINE_TODO 已插入迁移提示 -> "
          f"{len(io.open(P, encoding='utf-8').read().encode('utf-8'))} B")
else:
    print("⚠ 未找到「接手入口」标题")
