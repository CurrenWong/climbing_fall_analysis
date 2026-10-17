"""固化 CART_PATCH 一轮：4 条结论 + 2 个真 bug 修复。"""
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
L = ROOT / ".learnings/LEARNINGS.md"

REC = """
---

## LRN-20261004-053 [mixed] 软骨重分网打通「可建可解」，但 tie 界面成为新瓶颈

**任务**：按 `[LRN-052]` 的结论（细分/平滑/改方向三条捷径全部排除）**重新分网关节底面**。

### ✅ 成果：细网格已可建、可解（rc=0）
新增 `scripts/ankle_fe/cart_patch.py`（**2D 参数平面重心网格 + 精确抬回**）：
- 不做边界裁剪、不做点定位：在**每个原三角内部**打重心网格 `(i,j,k), i+j+k=m`
  ⇒ 域 = 原补丁，覆盖率天然 100%
- 用**同一组重心坐标**精确抬回 ⇒ 新点精确落在原（分片线性）曲面上 ⇒ **底面翘曲恒为 0**
- 位移场同样用重心坐标插值（`[LRN-052]` 验证过的正确做法）

| m | 底面三角 | 软骨 tet | 接触面(tib/tal) | 耗时 | rc |
|---|---|---|---|---|---|
| 基线 | 1119 | — | 152/200 | 18s | 0 |
| 1 | 1808 | 1074 | 152/200 | 59s | **0** |
| 2 | 5424 | 3222 | 456/600 | 90s | **0** |
| 3 | 10848 | 6444 | 915/1206 | 156s | **0** |
| 4 | 14464(=4²×809+...) | — | — | 218s | **0** |

**质量对比（归一化体积 nv=V/lmax³，正四面体 0.1179）**——细网格**不比基线差**：

| | nv 中位 | nv<1e-4 |
|---|---|---|
| 基线 t0d talus | 1.31e-02 | **1.29%** |
| m=1 | 9.86e-03 | 0.65% |
| m=2 | 1.75e-02（更好） | **0.65%**（基线一半） |
| m=3 | 1.79e-02（更好） | 0.60% |

### 🔴 三个真 bug（全部已修，且基线回归通过）

#### ① `build_feb.py:494` 的 `DET_MIN` 逻辑漏洞（**已修，基线受益**）
```python
if d > DET_MIN: keep
else: 翻转; if d2 > 0: keep   # ★ 正但小的单元被翻成负 ⇒ 当坏单元剔除 = 留洞
```
**正体积就是有效单元**（小 ≠ 无效）。修法：`if d > 0: keep`（只有 d≤0 才翻转）。
实测：**基线剔除 9 → 0**、「小但正」正好 9 个；`rc=0/17.1s` 不变 ✓。

#### ② 底面三角的对角线**必须与骨面一致**（这是最初那个「无法收敛」的真凶）
我为了优化最小内角，给 26/35 个 quad 换成 1-3 对角；但骨面用的是 **0-2 对角**
（`build_feb` 的 `HEX_T` 分解）⇒ **底面三角与骨面三角不重合 ⇒ tie 把节点拉错
⇒ 虚假刚度 ⇒ 单元翻转**。
实测：强制 `CP_DIAG=02` 后，接触面 **152/200**、tie primary **152/200**
与基线**完全一致**，`rc=0`。
⇒ **可复用判据：tied-elastic 的 primary 三角化必须与 secondary（骨面）的同源同对角。**

#### ③ 细网格必须用 `AUGLAG`，`PENALTY` 必失败（**决定性判据**）
| | m=2 PENALTY | m=2 AUGLAG |
|---|---|---|
| rc | 1 | **0** |
| 失败点 | t=0.0045 | — |

**判据**：失败点恒为 `1/(11·TIME_STEPS)`，与步长无关（TS=50/100/200 → 全是 /11）
⇒ **近乎零载荷就翻转** ⇒ 不是几何/载荷问题，**是接触求解器**。
另外 TIME_STEPS 扫描（50/100/200）**完全无法改善** ⇒ 排除「步长太大」。

### 🚧 剩余瓶颈：tie 界面
**现象**：
- `|u|max = 1.35~2.08 mm`，而指定的顶面位移只有 **0.1472 mm**（**10×以上**）
- 支反力**不收敛**：penalty=1 时 308.5→205.1→227.4 N（m=1/2/3）
- **penalty 1 能跑，4/16/64 全部 rc=1**（与 `[LRN-046~050]` 里 penalty 单调收敛的
  结论**相反** —— 那里是 `cart_refine` 的 L=0 粗网格）

⇒ 解释：`TIE_PENALTY=1`（手册默认）**太软**，软骨底面在滑移 ~1 mm
（`|u|max` 远超指定位移即是此特征），所以支反力没有意义；而加硬罚因子又发散。
**⇒ tie 界面（罚因子/搜索半径/secondary 覆盖）是细网格下的新瓶颈。**

### 📦 产物
- `scripts/ankle_fe/cart_patch.py`（新，`CP_M`/`CP_DIAG`/`CP_DOMAIN`/`CP_OUT`）
- `scripts/ankle_fe/build_feb.py`（修 `DET_MIN`，新增「小但正」计数）
- `temp/pyfebio_demo/`：`feas_patch2.py`（可行性）/`diag_degen.py`/`diag_nv_m.py`/
  `conv_final.py`/`pen_conv2.py`/`conv_study.py`
- 网格产物：`cartilage_mesh_patch_m{1,2,3,4}.npz`

### ➡️ 下一步
1. **查清 tie 为何在细网格下发散**：`TIE_SEC_TOL`（默认 0.05 mm）扩大？`search_radius`
   （1.0 mm）调整？secondary 覆盖是否够（m=2 时 primary 456 面 vs secondary 仅 117 面）？
2. 或**绕开 tie**：把骨关节面也细分到与软骨共形（此前「加密骨」失败，但那时是整根骨；
   只细分关节区局部可能可行）
3. 收敛判据应改用 **`|u|max ≈ 指定位移`** 作 sanity（现在 10× 明显异常）

- 日期：2026-10-04
"""
if "LRN-20261004-053" not in L.read_text(encoding="utf-8"):
    L.write_text(L.read_text(encoding="utf-8").rstrip() + "\n" + REC, encoding="utf-8")
print(f"LEARNINGS.md -> {L.stat().st_size:,} B")

P = ROOT / "CART_HEX_TODO.md"
P.write_text("""# CART_HEX — 关节底面重分网，提升接触分辨率

**目标**：关节面从 1119 个底面格 → 数千；并解开 P1d2b 的 hex8 翻负阻塞。

## 状态：🟡 **可建可解（rc=0）**，但 **tie 界面成为新瓶颈**（见 `[LRN-20261004-053]`）

### ✅ 已完成
- [x] 三条捷径（改偏移方向 / Taubin 平滑 / 纯细分）实测排除 → `[LRN-052]`
- [x] 可行性实测：补丁对最佳拟合平面**无折叠 / 单连通 / 浅度 0.10~0.16** → 平面参数化可用
- [x] **`scripts/ankle_fe/cart_patch.py`**：2D 参数平面重心网格 + 精确抬回
      （翘曲恒为 0；域 = 原补丁，覆盖 100%）
- [x] **修 `build_feb.py` 的 `DET_MIN` 逻辑漏洞** → 基线剔除 9→0，回归 rc=0/17.1s ✓
- [x] **找到对角线约束**：底面三角必须与骨面同对角（0-2），否则 tie 拉错
- [x] **找到 AUGLAG 必需**：细网格下 PENALTY 必失败（失败点恒 1/(11·TS)）
- [x] m=1/2/3/4 全部 rc=0；分辨率 4×/9×/16×；质量不比基线差

### 🚧 未完成（当前瓶颈）
- [ ] **tie 界面**：`TIE_PENALTY=1` 太软 ⇒ 软骨底面滑移 ~1 mm
      （`|u|max` 1.35~2.08 mm vs 指定 0.1472 mm，**10×以上**）
      ⇒ 支反力无意义（308/205/227 N 不收敛）；罚因子加硬到 4/16/64 又 rc=1
- [ ] 查 tie 发散原因：`TIE_SEC_TOL`(0.05) / `search_radius`(1.0) / secondary 覆盖
      （m=2 时 primary 456 面 vs secondary 仅 117 面）
- [ ] 备选：绕开 tie，把**关节区局部**骨面细分到与软骨共形

## 关键开关
| 变量 | 值 | 说明 |
|---|---|---|
| `CP_M` | 1..4 | 每三角边分段数（子三角 = m²） |
| `CP_DIAG` | **必须 `02`** | 底面三角对角必须与骨面一致 |
| `CP_DOMAIN` | `bary` | `delaunay` 实测更差（瘦长是补丁固有） |
| `LAUGON` | **必须 `AUGLAG`** | 细网格下 PENALTY 必失败 |
| `TIE` | `1` | 底面已复制，必须开 tie |
| `CART_MESH` | `cartilage_mesh_patch_m{m}.npz` | |

## 复现
```
set CP_M=2 & set CP_DIAG=02 & set CP_OUT=cartilage_mesh_patch_m2.npz
python scripts\\ankle_fe\\cart_patch.py
set CART_MESH=cartilage_mesh_patch_m2.npz & set TIE=1 & set LAUGON=AUGLAG
set DRIVER=disp & set DISP_MM=0.1472 & set TIME_STEPS=20
python scripts\\ankle_fe\\build_feb.py
```

## 验收判据
1. 关节面底格 ≥ 1000 ✓（m=2 → 5424）  2. 剔除 0 ✓  3. rc=0 ✓
4. **`|u|max` ≈ 指定位移**（0.1472 mm） ❌ 现在 1.35~2.08 mm
5. 应力/接触面积随 m **收敛** ❌ 现在不收敛
""", encoding="utf-8")
print(f"CART_HEX_TODO.md -> {P.stat().st_size:,} B")
