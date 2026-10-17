"""落盘本轮硬发现：ERR（载荷被丢）+ LRN（力控制失败真因 = 关节初始张开）。"""
import re
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
LT = ROOT / ".learnings/LEARNINGS.md"

t = LT.read_text(encoding="utf-8")
have = set(re.findall(r"\[LRN-(\d{8}-\d{3})\]", t))
NEW = next(f"20261005-{k:03d}" for k in range(58, 95) if f"20261005-{k:03d}" not in have)
print(f"⇒ LRN-{NEW}")

ENTRY = r"""

---

## [LRN-NEW] 踝韧带线：`DRIVER=force` 失败的真因 = **关节初始张开**（不是刚体模态）；且**载荷被丢**是真 bug

- **日期**：2026-10-05
- **类别**：`correction` / `knowledge_gap`
- **触发**：执行 `ANKLE_GOAL_TODO.md` Step 1~2（先验收、再改 BC）

### 🔴 发现 1（真 bug）：`ankle_lig_feb.py` 的 `<Loads>` **恒为空**

`scripts/ankle_fe/ankle_lig_feb.py` 第 865~873 行：

```python
    if not _use_force or not _faces or F <= 0:
        continue
        # traction = 力 / 面积 ...        ← 缩进 8，属 if 体内
        a = area_proj(...)                 ← 以下 7 行全为 **死代码**
        w(f'    <surface_load ...')        ← 永远不执行
```

`continue` 之后整段代码**缩进在 `if` 体内** ⇒ `DRIVER=force` 时**载荷被静默丢弃**。

**实测**（`TRUNC_MM=60 FIX_TALUS=transverse`，`MATSET=WK`）：
| | `<Loads>` 内容 |
|---|---|
| 修前 | **0 行**（空） |
| 修后 | 6 行，`traction=-7.821332,-1.294633,-5.402078` ×2 ⇒ \|tr\|=**9.594 N/mm²** |

> 两张 tractions 相同是**正确的**：`FIB_SHARE=-1` 按面积分配（胫骨 47.8 / 腓骨 14.7 mm²
> ⇒ 腓骨 23.5%），故单位面积牵引相等。

**修法**：把第 867~873 行 dedent 4 空格到 `for` 体内（`edit_file` 已应用；`ast.parse` OK）。

**⚠️ 由此：此前文档里的「`DRIVER=force` 病态 = 胫骨只靠接触约束 ⇒ 刚体模态」是
不完整诊断** —— 第一因是**载荷压根没施加上**。凡"某驱动方式不行"的结论，**先验产物
（`.feb`/载荷段/边界段）再看求解器报错**。

### 🔴 发现 2：修好载荷后，`DRIVER=force` **仍然失败**，但与载荷大小无关

修后实测（`FIX_TALUS=all`，排除距骨 BC 干扰）：

| 配置 | 结果 |
|---|---|
| `DRIVER=disp DISP_MM=0.20 FIX_TALUS=all` | **NORMAL TERMINATION**，FEBio 15 s ✅（对照，文档声称的可用配置，复现成功） |
| `DRIVER=force FORCE_N=600 FIX_TALUS=all` | **43 negative jacobians**，首步 t=0.0045，4 s ❌ |
| `DRIVER=force FORCE_N=0.6` | **不收敛**，t=0.1545，**463 s** ❌ |
| `DRIVER=force FORCE_N=600 FIX_TALUS=transverse` | 46 negative jacobians，3.4 s ❌ |

**⇒ `FIX_TALUS=all` 救不了 ⇒ 距骨 BC 不是原因；`FORCE_N=0.6`（1/1000）也失败
⇒ 与载荷大小无关 ⇒ 模型**本身**病态，不是"载荷太大压翻"。**

### 🔴 发现 2 的真因（机制）：**关节初始是张开的**

- `[LRN-20261005-056]` 实测：t=0 时间隙**中位 0.9732 mm**、接触面积仅 **31.0 mm²**
  （加载压 0.2 mm 后才到 ~300 mm²）
- `TOPFIX=lateral` 下胫骨顶面**轴向自由**（力控制必需）；接触又几乎不承载
  ⇒ **轴向近刚体模态 + 病态接触** ⇒ 小载荷也磨不动、大载荷直接翻单元
- `DRIVER=disp` 之所以"能跑通"，是因为**位移被强制**，绕过了"关节未闭合"这一步
  ⇒ **它的 rc=0 不代表力控制可行**

**⇒ 这正是 Anderson 2006 §2.2 的「弱弹簧 + 200 N 落位」要解的两件事**（机制完全对上）：
- **弱弹簧 500 N/mm** ⇒ 给张开阶段的胫/距骨一个**确定的轴向参考与载荷通路**
- **200 N 落位** ⇒ 先闭合 0.97 mm 间隙、把接触斑从 31 mm² 扩到 ~300 mm²，再上正式载荷

### 🔴 发现 3：韧带线用的是**旧软骨网格**

构建日志（`MATSET=WK TRUNC_MM=60`）：

```
!!  tibiotalar_tibia  tet4 保留 456 剔除 6  det(J) min=5.0995e-01
!!  tibiotalar_talus  tet4 保留 596 剔除 16 det(J) min=5.0800e-01
翻转修正 291 个；仍无法修正而剔除 22 个        ← 软骨网格自身有 291 个翻转单元
```

**而 `build_feb.py` 那条线早已用 `cart_patch.py` 修好**（`[LRN-051]`：2D 参数平面重心网格
+ 精确抬回，底面翘曲恒 0，m=1..4 全 rc=0）。
⇒ **两条线不是同一份软骨网格**；韧带线沿用了旧的、含 291 个翻转 + 22 个被剔除（= 留洞）的网格。
`CART_DIR` 环境变量可指向 `cart_patch.py` 的输出。

### 📦 下一步（顺序已由本轮实验确定）

1. ✅ 修载荷 bug（已完成，`[ERR-20261005-007]`）
2. 把韧带线的软骨换成 `cart_patch.py` 的网格（`CART_DIR`）—— 先清掉 291 翻转 / 22 洞
3. 加 `TAL_SPRING`（距骨平移弱弹簧 500 N/mm）+ `REL_IE_IR`（落位后解除内/外翻、内/外旋）
4. 加**落位阶段**（轻接触 → 200 N → 正式载荷）
5. 再用 `DRIVER=force` 验收敛

### 🧠 方法论教训

- **"某配置不行"必须附产物证据**：载荷段 / BC 段 / 初始间隙。本轮第一因（载荷被丢）
  藏在 `.feb` 里，光看求解器报错永远看不出来。
- **`rc` 要取对来源**：`ankle_lig_feb.py` 用 `subprocess.run` 跑 FEBio 但**不透传 returncode**
  ⇒ 外层 `rc=0` 与 FEBio 是否成功**无关**，必须解析日志里的 `NORMAL/ERROR TERMINATION`。
- **"收敛"≠"物理可行"**：`disp` 收敛是因为绕过了关节张开；判断能力要看着**驱动力式**是否为目标所需。
- **同项目两套脚本可能各带一份资产**（此处=软骨网格），改一条线要查另一条的资产是否更新。
""" .replace("LRN-NEW", f"LRN-{NEW}")

LT.write_text(t.rstrip() + ENTRY, encoding="utf-8")
print(f"LEARNINGS.md -> {LT.stat().st_size:,} B  [{'OK' if NEW in LT.read_text(encoding='utf-8') else 'XX'}]")

# ── TODO 更新 ──
TD = ROOT / "ANKLE_GOAL_TODO.md"
s = TD.read_text(encoding="utf-8")
s = s.replace(
    "- [ ] **Step 1** 写验收脚本 `scripts/ankle_fe/acceptance.py`；对**当前配置**跑出\"失败\" ⇒ 自证判据有区分度",
    "- [x] **Step 1a** 核验产物：`ankle_lig_feb.py` 的 `<Loads>` **恒为空**（死代码）⇒ **载荷被丢** ⇒ 已修 "
    f"`[ERR-20261005-007]`（修后 0→6 行，|tr|=9.594 N/mm²）\n"
    "- [x] **Step 1b** 对照实验：`disp` **NORMAL 15s** ✅ vs `force` 600N **43 负 Jacobian 4s** ❌ / "
    "`force` **0.6N 463s 不收敛** ❌ ⇒ **与载荷大小无关，真因=关节初始张开**（中位间隙 0.97mm、t=0 接触 31mm²）"
    f"`[LRN-{NEW}]`\n"
    "- [x] **Step 1c** 发现韧带线用**旧软骨网格**（291 翻转 / 22 剔除留洞）；`cart_patch.py` 的修好版未被引用\n"
    "- [ ] **Step 1d** 写验收脚本 `scripts/ankle_fe/acceptance.py`；对当前配置跑出\"失败\" ⇒ 自证判据有区分度\n"
    "- [ ] **Step 2a** 韧带线软骨换 `cart_patch.py` 网格（`CART_DIR`）—— 先清掉 291 翻转/22 洞\n"
    "- [ ] **Step 2b** 加落位阶段（轻接触 → 200 N → 正式载荷）")
s = s.replace("- [ ] Step 2 `ankle_lig_feb.py` 加 `TAL_SPRING`", "- [ ] Step 2c `ankle_lig_feb.py` 加 `TAL_SPRING`")
s = s.replace("## 当日记录\n- 2026-10-05 建本清单；开始 Step 1",
              f"## 当日记录\n- 2026-10-05 建本清单\n- 2026-10-05 Step 1 完成：**载荷被丢**（真 bug，已修）+ "
              f"**力控制失败真因=关节初始张开**（与载荷大小无关，0.6N 也失败）⇒ `[LRN-{NEW}]` / `[ERR-20261005-007]`")
TD.write_text(s, encoding="utf-8")
print(f"ANKLE_GOAL_TODO.md -> {TD.stat().st_size:,} B")
print(f"⇒ LRN-{NEW}")
