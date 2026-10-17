"""① 撤回 MEMORY 里被错标实验支撑的「44.8% 已否证」；② 记录真实 A1+A2 矩阵结果。"""
import io
import os

W = r"C:\Users\curre\.copaw\workspaces\default"
R = r"D:\Project\climbing_fall_analysis"

# ---- MEMORY ----
MP = os.path.join(W, "MEMORY.md")
m = io.open(MP, encoding="utf-8").read()
OLD = ("- 🔍 **2~4× 的来源**：① ~~缺距下关节 44.8%~~ **已否证**（2026-10-06 实测：接入后 F "
       "**+0.0000%**，因 `FIX_TALUS=transverse` 冻结距骨⇒接触不承力）② 无摩擦切向滑移 "
       "**39.5%**（只 0.029mm 真压扁）③ tie 罚滑移 **13.9%**；骨 <6%")
NEW = ("- 🔍 **2~4× 的来源（⚠️ 44.8% 那条未定）**：① **缺距下关节 44.8% —— 未证实也未否证**"
       "（2026-10-06 的「+0.0000%」测量**因 case 未显式传 `SUBTALAR`、实际没接上而作废**；"
       "接入后能否收敛**尚未解决**）② 无摩擦切向滑移 **39.5%** ③ tie 罚滑移 **13.9%**；骨 <6%")
if OLD in m:
    m = m.replace(OLD, NEW)
    print("✓ MEMORY 已替换 44.8% 那条")
else:
    print("⚠ MEMORY 未命中（可能已被改），输出现值以便核对：")
    for ln in m.split("\n"):
        if "44.8" in ln:
            print("   ", ln[:160])

OLD2 = ("- **⏭️ 下一步**：① **A1+A2 一起做**（接入距下关节 **且** 解除 `FIX_TALUS=transverse`）"
        "—— 单接关节 **零效应**（实测 +0.0000%），因横向硬约束冻结距骨、接触不承力；"
        "而**接入后收敛需 `MAX_UPS=0`**（根因=拟牛顿更新病态，非几何）② 力控找根接 OpenSim")
NEW2 = ("- **⏭️ 下一步**：① **A1+A2 已跑 4 格矩阵**（E0/E1/E2/E3）⇒ 基线 592.6020 ✓，"
        "E1/E2/E3 **均无负 J、但在 t=0.475/0.585/0.95 收敛停滞（残差卡 ~1e2）**；"
        "`MAX_UPS=0` 已消掉 t=0 失控 ⇒ **下一步是迭代/步长**（`TIME_STEPS`↑、`MAX_REFS`↑）"
        "② 力控找根接 OpenSim")
if OLD2 in m:
    m = m.replace(OLD2, NEW2)
    print("✓ MEMORY 下一步已替换")
else:
    print("⚠ MEMORY 下一步未命中")
io.open(MP, "w", encoding="utf-8").write(m)
n = len(m.encode("utf-8"))
print(f"MEMORY {n} B " + ("✅" if n <= 10240 else f"❌超 {n-10240}"))

# ---- ROADMAP 追加更正 ----
ADD = """
---

## 🚨 更正（2026-10-06 深夜）：上一节「A1 真实结果」有两条建立在**错标实验**上

**根因**：`a1_matrix2.py` 的 case `"both_noupd": {"MAX_UPS": "0"}` **没有 `SUBTALAR` 键**，
而我是在提交 `85b72f1`（把默认从 `auto` 改成 `0`）**之后**跑它的
⇒ 那个叫「双侧全接」的 case 实际是 **`SUBTALAR=0`、根本没接上**。
`a1_isolate.py` 的 `sub_noupd_ctl` 同样没传 ⇒ 也未接上。

**必须撤回的两条**：
1. ~~「接入距下关节的真实效应 = +0.0000%」「44.8% 被直接测量否证」~~
   ⇒ **作废**。真实数字**未知**（接入后尚未收敛到 t=1）。
2. ~~「`MAX_UPS=0` 后双侧接入完全收敛 rc=0」~~ ⇒ **作废**（那格没接上）。

**仍然成立的**：
- `MAX_UPS=0` **确实消掉了 t=0 的失控**：接入配置从「11 负 J、t=0 步」变为 **0 负 J、走到 t=0.475**
  ⇒ H1（拟牛顿更新病态）**方向正确，但不充分**。
- 几何预检结论（执行者交付）不变。

### ✅ 真实 A1+A2 矩阵（`temp/pyfebio_demo/a1a2_matrix.py`，显式传参）

| case | 配置 | rc | 负J | 失收敛 t | 方程数 | 两顶面合力 |
|---|---|---|---|---|---|---|
| **E0** | 基线（不接 + fix） | 0 | 0 | —（收敛） | 11,112 | **592.6020 N** |
| **E1** | A1+`MAX_UPS=0` | 1 | **0** | **0.475** | 13,830（+2718=距下软骨） | — |
| **E2** | A1+A2+`MAX_UPS=0` | 1 | **0** | **0.585** | 14,476（再+646=接地弹簧） | — |
| **E3** | 仅 A2+`MAX_UPS=0` | 1 | **0** | **0.95** | 11,758（+646） | — |

**判读**：
- **接入确实生效**（方程数 +2,718 / +646 精确对应软骨与弹簧数量）
- **失败已从「发散/负 J」降级为「收敛停滞」**（残差卡在 ~1e2，基线收敛到 1.42e2）
- ⇒ **不是病态了，是迭代次数/时间步不够** ⇒ 下一步试 `TIME_STEPS`↑（20→40）
  与 `MAX_REFS`↑（50→100），**而非继续怀疑几何**
"""
RP = os.path.join(R, "ANKLE_ROADMAP_TODO.md")
r = io.open(RP, encoding="utf-8").read()
if "🚨 更正（2026-10-06 深夜）" not in r:
    io.open(RP, "a", encoding="utf-8").write(ADD)
    print(f"✓ ROADMAP 追加 {len(ADD)} 字符 -> "
          f"{len(io.open(RP, encoding='utf-8').read().encode('utf-8')):,} B")
