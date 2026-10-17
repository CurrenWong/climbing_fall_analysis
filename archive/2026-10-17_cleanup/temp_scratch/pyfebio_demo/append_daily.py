"""把本轮 TIE 成果追加到 daily note。"""
import io
import os

DST = r"C:\Users\curre\.copaw\workspaces\default\memory\2026-10-05.md"
SEP = "\n\n---\n\n"
BLOCK = """## ✅ 踝子线：TIE 移植成功 ——「~235 N 天花板」被打破（`[LRN-20261005-064]`）

**用户要求**：① 给 `ankle_lig_feb.py` 移植 `TIE` ⇒ 用 `cart_patch.py` 的干净软骨网格（治 ~235 N 上限）；
② 力控制（Anderson 弱弹簧）⇒ 接上 OpenSim。**本轮完成 ①，② 未开始。**

### 硬结果（`TIE=1 CART_NEW=1 CART_DIR=...cart_m1 DET_MIN=0 ALIGN_AXIS=1`）

| DISP (mm) | 0.25 | 0.50 | 1.00 | 1.50 | **1.60** | 1.70 |
|---|---|---|---|---|---|---|
| Fz (N) | 124.0 | 230.9 | 418.1 | 567.2 | **592.6** | ERROR |
| 状态 | NORMAL | NORMAL | NORMAL | NORMAL | **NORMAL** | ✗ |

**位移窗口 0.3 → 1.5~1.6 mm（5×）；力 ≤235 → 592.6 N（2.5×）。600 N 达标。**
旧格式**零回归**：`209.5589 N`(0.25) / `183.2044 N`(0.20) 逐字一致。

### 不确定性带（`[LRN-062]` 自查，必须随数字报）
- 关节接触罚因子**不是**柔度来源（`PENALTY 0.1→1` 仅 −8.5%；AUGLAG +4.2%）✅
- **tie 罚因子是**：PEN=1/3/10/30 → 400.8/420.2/449.1/494.3 N ⇒ **带 ~+7~18%（记 ±20%）**

### 修的 4 个真 bug
1. `{j}_elems`/`{j}_tie_faces` 用**负索引** `-k` ⇒ 行 `k-1`（我按 `x-1` 写 ⇒ `KeyError`）
2. tie 接触**必须写 `<laugon>PENALTY</laugon>`**（漏掉 ⇒ t=0 负 J，与位移无关）
3. tie 的 `<contact>` 要写在关节接触**之前**
4. ⭐ **primary 反查必须同时扫 `tet4` 和 `penta6`** —— 我改成 penta6 后反查静默取空 ⇒
   **tie 被跳过、软骨脱开**，日志仅一句"跳过" ⇒ 我据此误判数轮

### 方法论教训
- **`acceptance.py` 的「负 Jacobian N 处」是假信号**：`len(findall("negative jacobian"))` = 消息**行数**，
  非单元数 ⇒ 必须读 FEBio 原文 `* NN negative jacobians detected. *`
- **cmd `set B=...& %B%` 在 `for` 里不可靠** ⇒ 静默跑成默认配置，造出一堆假失败 ⇒ 改用 Python 驱动
- **剔单元往往比留坏单元更糟**：`DET_MIN=0.5` 剔 153 个 ⇒ 46 个负 J；`DET_MIN=0` ⇒ 只剩 2 个
- **多点采样才查得出 penta6 负 J**（质心采样会因正负抵消漏判），但顶点处 det 可合法逼近 0 ⇒
  **判据以 FEBio 输出为准，自检只用于定位**
- **m=1 收敛，m=2/3/4 在 t=0 报负 J** ⇒ 加密受限，须先修 `cart_patch.py` 的无效楔形

### 新增/改动的产物
- `scripts/ankle_fe/ankle_lig_feb.py`：新格式软骨加载（`CART_NEW`/`CART_MESH`/`CART_DIR`）
  + `TIE=1`（primary 反查存活单元、secondary KD-tree 就近骨面、`<tied-elastic>` 写入）
- `temp/pyfebio_demo/`：`curve_m1.py`（力-位移曲线）、`pen_check_m1.py` / `tie_pen_scan.py`（罚因子检验）、
  `bracket600.py`、`chk_penta_J.py` / `prism_thick.py` / `chk_twist.py`（网格自检）、`run_m.py`
- `.learnings/LEARNINGS.md` +LRN-064（151,883 → 156,833 B）
- `ANKLE_FIX_TODO.md` / `MEMORY.md` / 本条 daily

### 下一步
**力控制**：600 N 只需 ~1.62 mm（在稳定窗内）⇒ 可用「位移驱动 + 找根」，
不一定非要 Anderson 的 `<Discrete>` 弱弹簧。
"""
if "[LRN-20261005-064]" in io.open(DST, encoding="utf-8").read():
    print("⚠ 本节已存在，跳过")
else:
    b = os.path.getsize(DST)
    io.open(DST, "a", encoding="utf-8").write(SEP + BLOCK)
    print(f"daily {b} -> {os.path.getsize(DST)} bytes")
