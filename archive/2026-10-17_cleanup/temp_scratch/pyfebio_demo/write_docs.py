"""收口：更新 CART_HEX_TODO / ANKLE_FIX_TODO / 总体进度 / MEMORY / TIE_DIVERGE_TODO。"""
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")

# ── ① CART_HEX_TODO.md 重写 ──
(ROOT / "CART_HEX_TODO.md").write_text("""# CART_HEX — 关节底面重分网

## 🚨 状态修正（2026-10-05，见 `[LRN-20261005-054]`）

**「tie 是瓶颈」是误诊，已撤回。** 真根因在**关节接触的几何**：

- **关节两软骨面不共形** ⇒ 中位间隙 **1.6 mm**、只有 **1~3 个节点**接触（<0.05 mm）
  ⇒ **关节从未闭合** ⇒ **F 完全是罚因子伪影**：
  | 关节 PENALTY | 1e-4 | 1e-3 | 1e-2 | 1e-1 | 1 |
  |---|---|---|---|---|---|
  | PENALTY 模式 F (N) | 0.83 | 6.99 | 50.65 | 220.08 | **419.56** |
  | AUGLAG 模式 F (N) | — | 2.00 | 15.31 | 110.63 | **315.06** |
- ⇒ `436.8 N` / `419.56` / `190.80` / `315.06` **都不是收敛解**；
  「外推 600 N 胫侧 ≈3.4 MPa 对吴恺 3.74 MPa 吻合」**失效**
- **`|u|max` 不是判据**（胫骨 400 mm 悬臂，刚体拟合 0.17° 占位移 79%）
- **tie 本身没问题**：界面滑移中位 0.048 mm = 自身运动的 **3%**

## ✅ 已完成
- [x] 三条捷径（改偏移方向 / Taubin 平滑 / 纯细分）实测排除 → `[LRN-052]`
- [x] **`scripts/ankle_fe/cart_patch.py`**：2D 参数平面重心网格 + 精确抬回（翘曲恒 0）
- [x] 修 `build_feb.py` `DET_MIN` 逻辑漏洞 → 基线剔除 9→0
- [x] 找到 **CP_DIAG=02 铁律**（底面三角对角须与骨面同源）
- [x] 找到 **细网格须 `LAUGON=AUGLAG`**
- [x] m=1/2/3/4 全部 rc=0；质量不比基线差

## 🚧 未完成（真瓶颈）
- [ ] **修关节几何让两软骨面共形**（唯一真解）
      ① 统一用一张关节面 ② 或给初始过盈 0.05~0.1 mm ③ 或按 Anderson：骨刚性 + 只取关节局部
- [ ] 几何修好后**重跑罚因子收敛**（×10 与 ×0.1，F 变化 <2% 才算通过）
- [ ] 之后再谈 m2 vs baseline 的 1.65× 差异（现在几何都不同，无法归因）

## 关键开关
| 变量 | 值 | 说明 |
|---|---|---|
| `CP_M` | 1..4 | 每三角边分段数（子三角 = m²） |
| `CP_DIAG` | **必须 `02`** | 底面三角对角必须与骨面一致 |
| `CP_DOMAIN` | `bary` | `delaunay` 实测更差 |
| `LAUGON` | **必须 `AUGLAG`** | 细网格下 PENALTY 必失败（因关节罚太大） |
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

## 验收判据（已替换）
1. 关节面底格 ≥ 1000 ✓  2. 剔除 0 ✓  3. rc=0 ✓
4. ~~`|u|max ≈ 指定位移`~~ **作废**（长骨悬臂弯曲会合法超过）
5. **新**：`F` 对关节罚因子 ×10/×0.1 不敏感（<2%）
6. **新**：末态关节面**中位间隙 ≈ 0**、接触面 ≥ 30%（现在 **1.4%**）
""", encoding="utf-8")
print(f"CART_HEX_TODO.md -> {(ROOT/'CART_HEX_TODO.md').stat().st_size:,} B")

# ── ② TIE_DIVERGE_TODO.md 结案 ──
(ROOT / "TIE_DIVERGE_TODO.md").write_text("""# TIE_DIVERGE — 查 tie 界面发散的真因

## ✅ 已结案（2026-10-05）→ `[LRN-20261005-054]`

**结论：「tie 发散」是误诊。** tie 没有问题；真根因是**关节两软骨面不共形 ⇒ 关节从未闭合
⇒ F 完全是罚因子伪影**。

- [x] 1. 取地面真值（读 .feb 的 tie 块 / surface / NodeSet / BC）
- [x] 2. 复算覆盖：primary 456 面/313 节点、secondary 117 面/131 节点（能覆盖）
- [x] 3. 测实际间隙：**两软骨外面中位 1.6~2.4 mm，仅 1~3 点接触** ⇒ 几何根因
- [x] 4. `|u|max` 口径查清 = **全模型节点位移模最大值**；胫骨刚体拟合 **0.17° 转动占 79%**
      ⇒ 是 400 mm 悬臂弯曲，**不是 tie 滑移** ⇒ 旧判据作废
- [x] 5. 单变量扫：`TIE_SEC_TOL`/`search_radius`/tie `auto_penalty`/tie `penalty`
      **全部无关**；**关节 `PENALTY` 是唯一控制量**
- [x] 6. 结论 + 修复方向，实测验证（4 因子矩阵 + 罚扫描 + AUGLAG 扫描 + 位移扫描 + 载荷截面扫描）

### 关键实验（全部实跑）
| 实验 | 结果 |
|---|---|
| 4 因子矩阵（网格×TIE×LAUGON） | `LAUGON` 单独改变 F 25%；m2+PENALTY rc=1 |
| D 失败归因 | `Negative jacobian` @第一步；关 tie `auto_penalty`/加硬 tie penalty **无效** |
| 关节 `PENALTY=0.001` | ✅ rc=0 且自动步长主动放大 |
| 关节罚扫描 | F：0.83→6.99→50.65→220.08→**419.56 N**（1e-4→1） |
| AUGLAG 初始罚扫描 | F：2.00→15.31→110.63→**315.06 N**（1e-3→1）**也不收敛** |
| 位移扫描（P=0.001） | 2.0 mm 后**饱和在 ~80 N**，中位间隙停在 0.6 mm |
| 载荷截面 H=60/30 mm | 罚敏感比仍 67×/69×（`|u|max` 降了但 F 没救） |

### 下一步
1. **修关节几何让两软骨面共形**（唯一真解）
2. 重建罚收敛：`F` 对 ×10/×0.1 不敏感才算通过
3. 之后重做 m2 vs baseline 对比
""", encoding="utf-8")
print(f"TIE_DIVERGE_TODO.md -> {(ROOT/'TIE_DIVERGE_TODO.md').stat().st_size:,} B")

# ── ③ ANKLE_FIX_TODO.md 顶部插警示 ──
AF = ROOT / "ANKLE_FIX_TODO.md"
s = AF.read_text(encoding="utf-8")
BAN = """
> ## 🚨 2026-10-05 重大更正（见 `[LRN-20261005-054]`）
> **本文件里所有 `F` 数值都不是收敛解。** 实测：关节两软骨面**不共形**（中位间隙 1.6 mm、
> 仅 1~3 点接触）⇒ **关节从未闭合** ⇒ `F` 完全由关节接触罚因子决定：
> `PENALTY` 1e-4→1 时 **0.83 → 419.56 N**，AUGLAG 也不收敛（2.00 → 315.06 N）。
> ⇒ **`436.8 N` 与「外推 600 N 胫侧 ≈3.4 MPa 对吴恺 3.74 MPa 吻合」失效。**
> ⇒ 唯一真解：**修关节几何让两软骨面共形**，再重跑罚收敛。
"""
if "2026-10-05 重大更正" not in s:
    lines = s.split("\n")
    AF.write_text(lines[0] + "\n" + BAN + "\n".join(lines[1:]),
                  encoding="utf-8")
print(f"ANKLE_FIX_TODO.md -> {AF.stat().st_size:,} B  "
      f"[{'OK' if '2026-10-05 重大更正' in AF.read_text(encoding='utf-8') else 'XX'}]")

# ── ④ docs/总体进度.md 子线表格加警示行 ──
OP = ROOT / "docs/总体进度.md"
s = OP.read_text(encoding="utf-8")
OLD = "| **新** 关节底面重分网 `cart_patch.py` | 🟡 |"
NEW = ("| **🚨 关节几何不共形** | 🔴 | **F 全是罚因子伪影**（PENALTY 1e-4→1 ⇒ 0.83→419.56 N；"
       "AUGLAG 也不收敛 2.00→315.06 N）⇒ 关节从未闭合（中位间隙 1.6 mm、1~3 点接触）"
       "⇒ `436.8N`/「3.4≈吴恺3.74」**失效**。`[LRN-054]` |\n"
       "| **新** 关节底面重分网 `cart_patch.py` | 🟡 |")
if OLD in s and "LRN-054" not in s:
    OP.write_text(s.replace(OLD, NEW, 1), encoding="utf-8")
print(f"docs/总体进度.md -> {OP.stat().st_size:,} B  "
      f"[{'OK' if 'LRN-054' in OP.read_text(encoding='utf-8') else 'XX'}]")

# ── ⑤ MEMORY.md 脚踝段更正（等量/更小） ──
M = WS / "MEMORY.md"
s = M.read_text(encoding="utf-8")
OLD = """- ✅ tied-elastic（`TIE=1`）L=0 单调收敛 436.817→432.5 N（−0.99%）；
  🟡 **软骨已可重分网**（`cart_patch.py`，m=1..4 全 rc=0，底面格 5424@m2）
  ⚠️ 铁律：**底面三角对角必须=骨面(0-2)**，否则 tie 拉错⇒假刚度⇒翻转；
  **细网格必须 `LAUGON=AUGLAG`**。🚧 瓶颈=tie 界面（`|u|max` 超指定 10×、F 不收敛）
  `[LRN-052/053]`"""
NEW = """- 🚨 **所有 F 都不是收敛解！** 关节两软骨面**不共形**（中位间隙 1.6mm、仅 1~3 点接触）
  ⇒ 关节从未闭合 ⇒ F 全由关节罚因子定：PENALTY 1e-4→1 得 **0.83→419.56 N**，
  AUGLAG 也不收敛（2.00→315.06 N）⇒ `436.8N`/「3.4≈吴恺3.74」**失效**。
  真解=**修几何让软骨面共形**再重跑罚收敛。`|u|max` **非判据**（胫骨 400mm 悬臂弯 0.17°）
- ✅ tie **没问题**（界面滑移=自身运动 **3%**）；软骨可重分网（`cart_patch.py`，m=1..4 rc=0）。
  铁律：**底面三角对角须=骨面(0-2)**；**细网格须 `LAUGON=AUGLAG`**。`[LRN-052/053/054]`"""
if OLD in s:
    M.write_text(s.replace(OLD, NEW, 1), encoding="utf-8")
sz = M.stat().st_size
print(f"MEMORY.md -> {sz:,} / 10240  余量 {10240 - sz:+,}  "
      f"[{'OK' if 'LRN-054' in M.read_text(encoding='utf-8') else 'XX'}]")
