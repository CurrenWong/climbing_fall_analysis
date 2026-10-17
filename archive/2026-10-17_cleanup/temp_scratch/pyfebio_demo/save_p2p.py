"""落盘 P2' 发现：加载面积不一致(ERR) + 韧带 EA 纠正 + 计划更新。"""
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"
TOD = os.path.join(R, "ANKLE_FIX_TODO.md")
ERR = os.path.join(R, ".learnings", "ERRORS.md")
DAY = os.path.join(W, "memory", "2026-10-06.md")

# ── ERR 条目（扫号防撞）──
et = io.open(ERR, encoding="utf-8").read()
hdr = re.findall(r"^#+\s*\[?(ERR-\d{8}-\d{3})", et, re.M)
tod_err = sorted({int(h.split("-")[2]) for h in hdr if "-20261006-" in h})
print(f"今日 ERR 标题号: {tod_err or '无'}")
nxt = (max(tod_err) + 1) if tod_err else 1
print(f"⇒ 下一个 ERR 号: ERR-20261006-{nxt:03d}")

EADD = f"""

---

## [ERR-20261006-{nxt:03d}] 加载面积 `top_area` 与几何不一致（11~45×）⇒ 力控模式结果不可信

**日期**: 2026-10-06
**项目**: climbing_fall_analysis / 踝关节 FE
**严重度**: 高（影响一类实验的全部结论）

### 现象
`ankle_lig_feb.py` 打印「顶面投影面积：胫骨 47.8 mm²，腓骨 12.7 mm²」，
`build_feb.py:552` 用 `traction = FORCE_N / top_area` 施加力控载荷。

### 实测几何（新脚本 `temp/pyfebio_demo/p2p_toparea.py`）
| Surface | 三角 | 真实面积 | 投影到加载轴(z) | max n·z 均值 |
|---|---|---|---|---|
| `tibia_top` | 114 | **2148.3 mm²** | **555.5 mm²** | 0.349 |
| `fibula_top` | 33 | **745.7 mm²** | **108.3 mm²** | 0.343 |

⇒ 模型用的 47.8 mm² 与投影面积差 **11.6×**、与真实面积差 **45×**。

### 根因（待确认）
`tibia_top` 的 **z 范围 = −669.85 … −640.28（跨 29.6 mm）** ⇒ 它不是切平面，
而是一条**外表面带**（含内踝区）。`area_proj()` 可能只累加了近轴面片
（|n·z| 阈值），因此得到 47.8 而不是投影面积 555.5。
**⇒ 待办：读 `area_proj` 与 `top_tibia` 的选取逻辑，确认哪个面积才是 BC 实际作用面。**

### 影响（关键，必须区分两类实验）
- **❌ 力控（`DRIVER=force`）**：除数错 11~45× ⇒ **此前所有力控实验结果不可信**
  （含"斜坡力控制"、`SEAT_HOLD`、200 N 预载等）。**必须在修好面积后重做。**
- **✅ 位移控制（`DRIVER=disp/dispx`）**：BC 是 prescribe 到 `tibia_top_ns`
  的 165 个节点，**位移值精确** ⇒ **P1 骨折裕度表、C2/C3、压缩曲线不受影响。**

### 修复方向
1. 修 `area_proj` 或改用**真实承载截面**（切面 ≈152~264 mm²，据 `p2p_crosssec.py`）
2. 加一条自检：**打印的面积必须与几何凸包/切面面积一致（容差 10%）**，
   不一致就报错而不是静默继续
3. 修完后重跑力控实验
"""
if f"ERR-20261006-{nxt:03d}" in et:
    print("· ERR 已有")
else:
    io.open(ERR, "a", encoding="utf-8").write(EADD)
    print(f"ERRORS +{len(EADD)} 字符 -> "
          f"{len(io.open(ERR, encoding='utf-8').read().encode('utf-8'))} B")

ADD = """

### 🔬 2026-10-06 第九轮：**P2' 韧带刚度（进行中）—— 抓到一个真 bug 级不一致**

#### ① 纠正我自己：壳韧带 EA 在解剖范围内（上一个判断读错了列）
`lig_read_*.npz` 里的 `lig_ea`（180/510/75.4）是 `ligaments.py` 的 **`EA(MPa)` 元数据列**，
**shell 模式不用它** —— `LIG_E_UNIFORM=True`（MATSET=WK）时壳统一 `E=260`，
真实 **EA = 260 × 截面**（`ligaments.py` 的 `A_cross = 膜面积/L × 厚`）：

| 韧带 | 截面(mm²) | 壳 EA (N) | 解剖 EA | 判定 |
|---|---|---|---|---|
| ATFL | 29.05 | 7,553 | 500~3600 | 偏硬 |
| CFL | 17.10 | 4,446 | 1000~9000 | ✓ |
| PTFL | 14.84 | 3,858 | 1000~7500 | ✓ |
| MD1/2/3/4 | 17~24 | 4,155~6,175 | 1000~12000 | ✓ |
| PTCL/SYND/SPRING | 25.9~62.0 | 6,739~16,130 | — | — |

**⇒ 之前"偏软 3~40×"的结论撤回**（那是把元数据列当成了实际刚度）。

#### ② 🚨 抓到 `top_area` 与几何差 11~45×（→ `[ERR-20261006-003]`）
见 ERRORS.md。**⇒ 所有力控实验结果作废待重做；位移控制结果不受影响** ✓

#### ③ 仍未解：τ 柔度归属（`370 N/mm` 的出处）
当前账对不上：若关节 k≈42,000 且韧带 k_axial≈285 ⇒ 合计 ≈10,000 ≫ 实测 370。
⇒ 差一个 ~384 N/mm 的**串联**柔度。**头号嫌疑已转到 `tibia_top` 那条 29.6 mm
高的外表面带**（BC 作用面不是切面）⇒ 局部弯曲/载荷扩散。
**⇒ 下一步（P2'' ）**：修/确认 `area_proj` 与 BC 作用面，再重测刚度。

#### 计划更新
- **P2' 拆成两件**：**P2'a 修加载面积 + 自检**（先做，它挡住力控）｜
  **P2'b 柔度归属**（`tibia_top` 带 vs 切面）。
- **P4 力控找根被 P2'a 阻塞**（面积不修，力控无意义）。
- **P1 骨折裕度表保持有效** ✓（位移控制）。
"""
t = io.open(TOD, encoding="utf-8").read()
if "第九轮：**P2' 韧带刚度" in t:
    print("· TODO 已有")
else:
    io.open(TOD, "a", encoding="utf-8").write(ADD)
    print(f"ANKLE_FIX_TODO +{len(ADD)} -> "
          f"{len(io.open(TOD, encoding='utf-8').read().encode('utf-8'))} B")

D = """
## 🔬 P2' 韧带刚度（第九轮）：纠正一条 + 抓到一个真不一致

**① 纠正**：`lig_ea`(180/510/75.4) 是元数据列，**shell 模式不用**；真实 EA=260×截面
⇒ ATFL 7,553 / CFL 4,446 / PTFL 3,858 N，**在解剖范围内** ⇒ 我上一条"偏软 3~40×"**撤回**。

**② 🚨 抓到 `top_area`**：模型用 47.8 mm²，实测投影 **555.5**、真实 **2148**
（`tibia_top` 是跨 29.6mm 的**外表面带**，不是切面）⇒ 差 **11~45×**
⇒ **力控结果全部作废待重做**；**位移控制不受影响（P1 骨折表有效）**。

**③ 未解**：370 N/mm 的柔度归属 —— 关节 42,000 + 韧带 285 ≈ 10,000 ≫ 370
⇒ 缺一个 ~384 N/mm 串联柔度，嫌疑在 `tibia_top` 那条带（BC 非切面）。
"""
d = io.open(DAY, encoding="utf-8").read()
if "P2' 韧带刚度（第九轮）" in d:
    print("· 日记已有")
else:
    io.open(DAY, "a", encoding="utf-8").write(D)
    print(f"日记 +{len(D)} 字符")
