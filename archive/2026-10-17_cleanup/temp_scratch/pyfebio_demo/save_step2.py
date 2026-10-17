"""落盘 Step 2 成果：稳定姿势配置 + 读数出口 + 验收脚本。"""
import re
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
LT = ROOT / ".learnings/LEARNINGS.md"

t = LT.read_text(encoding="utf-8")
have = set(re.findall(r"\[LRN-(\d{8}-\d{3})\]", t))
NEW = next(f"20261005-{k:03d}" for k in range(59, 99) if f"20261005-{k:03d}" not in have)
print(f"⇒ LRN-{NEW}")

ENTRY = r"""

---

## [LRN-NEW] 踝子线：**`DRIVER=dispx` + `FIX_TALUS=transverse` 是稳定的姿势可用配置**；并补上文本读数出口与验收脚本

- **日期**：2026-10-05
- **类别**：`best_practice` / `knowledge_gap`
- **承接**：`[LRN-20261005-058]`（力控制失败真因 = 关节初始张开；载荷被丢是真 bug）

### ✅ 结论一：稳定的**姿势可用**配置找到了

| 配置 | 结果 | FEBio 用时 | 顶面反力 | \|u\|max | 负 Jacobian |
|---|---|---|---|---|---|
| `DRIVER=disp`（顶面 xyz 全约束） | NORMAL ✓ | 15 s | — | — | 0 |
| **`DRIVER=dispx` + `FIX_TALUS=transverse`** | **NORMAL ✓** | **5 s** | **[90.0, 0, 0] N** | **0.382 mm** | **0** |
| `DRIVER=dispx` + `FIX_TALUS=all` | NORMAL ✓ | 31.7 s | [0, 0, 161.8] N | 0.535 mm | 0 |
| `DRIVER=force`（任意载荷，含 0.6 N） | **ERROR ✗** | 4~463 s | — | — | 41~46 |

**⇒ 三个要点：**
1. **`dispx + transverse` 两个都放开（顶面横向自由 + 距骨轴向自由）却依然收敛** —— 因为
   **位移控制不可能过冲**（这正是 `force` 失败而 `disp` 可用的同一机制，见 `[LRN-058]`）。
2. **`transverse` 不是 no-op**：与 `all` 相比反力 90 vs 162 N、`|u|max` 0.382 vs 0.535 mm，
   **响应显著不同** ⇒ 放松距骨 BC 确实改变了模型行为。
3. **`dispx` 有个必须记住的语义**：它 prescribe 的是 `kax = argmax(|LAX|)` 指定的**全局单轴**
   （本例 `LAX=[0.8153, 0.135, 0.5631]` ⇒ 只 prescribe **全局 x**），**不是真正的沿轴位移**。
   ⇒ 这解释了 `[LRN-038]` 里"`dispx` 是 no-op"的印象；用在姿势扫描前要先修正这一点。

### ⚠️ 结论二：一个**未能解释**的异常（必须留着，别当成已知）

`transverse` 的反力在 **x**，`all` 的反力在 **z** —— 但两者的顶面 BC 都是 prescribe 全局 x。
力平衡上，prescribe x 的约束反力**应有 x 分量**。**=> 待查**：
优先怀疑 ① 两次运行的 `reac_tibia` 节点集在 `.feb` 里是否真的一致；
② FEBio `Rx;Ry;Rz` 是否为"约束反力"还是"节点内力"；③ `disp` 段读数是否被上一步覆盖。
**在解释清楚前，不要用这两个数字做任何结算。**

### ✅ 结论三：补上了**文本读数出口**（此前只有二进制 xplt，无读者）

`scripts/ankle_fe/ankle_lig_feb.py` 新增 env 门控 `LOG_DATA`（默认 0 = 不改行为）。
语法出自 **FEBio 手册 §3.19.1 + 附录 E.2.1**（**不是猜的**）：

```python
LOG_DATA = os.environ.get("LOG_DATA", "0") == "1"
# <Output><logfile file="...">
#   <node_data data="ux;uy;uz" file="disp_<tag>.txt" delim=","/>
#   <node_data data="Rx;Ry;Rz" file="reac_<nm>_<tag>.txt" delim=",">节点列表</node_data>
```

**实测产出**（`WORK=temp/thums/febio_lig`，FEBio 4.13.0）：
`log_WK.txt`、`disp_WK.txt`（每步 3884 行 `id,ux,uy,uz`）、
`reac_tibia_WK.txt`（139 行）、`reac_fibula_WK.txt`、`reac_talus_WK.txt`（263 行）。
记录格式：`*Step = n` / `*Time = t` / `*Data = name` 后跟多列数据。

> **教训**：`ankle_lig_feb.py` 一直写着 `lig_read_{tag}.npz`「供 `read_lig.py` 用」，
> 但 **`read_lig.py` 从未存在**；而 `.xplt` 是**二进制**、无读者 ⇒ **历史上所有"读数"
> 都是临时脚本现算的，且已随 temp 清理丢失**。**这是"结论不可复现"的结构性原因** ——
> 现在有了正式出口。

### ✅ 结论四：验收脚本落地 → **Step 1d 完成**

`scripts/ankle_fe/acceptance.py`（新增，7.1 KB）：一条命令完成「跑 → 解析 → 出指标 → 判 C1」。
- 解析 `.feb`（节点坐标 / NodeSet / `<Loads>` 行数 / contact 设置）
- 解析 FEBio 文本读数（步记录）
- **Kabsch 刚体拟合**给转动角 + 残差占比；并单列"距骨低端刚体转动"
- 内置 `<Loads>` 空段回归闸门（`[ERR-20261005-007]` 的守门人）

用法：
```
python scripts/ankle_fe/acceptance.py DRIVER=dispx DISP_MM=0.20 FIX_TALUS=transverse
python scripts/ankle_fe/acceptance.py --no-run      # 只解析上一次结果
```

### 📌 与计划（`ANKLE_GOAL_TODO.md`）的差异

原计划 Step 2 是"加弱弹簧 `TAL_SPRING`"。**本轮实验后修正**：
- **`driving` 方式的选择比弹簧更根本**：换成 `dispx` 已经拿到稳定 + 放开转动，
  **不需要**先实现弹簧（弹簧是为了解决 `force` 的病态；而 `force` 现在可以不作为必经路径）
- 前提是接受「**用位移控制 + 对 DISP 做二分来匹配目标力**」（项目已有二分法传统）
  ⇒ 这样 **OpenSim 耦合仍然可达**，且绕开整个力控制病态问题
- **真正的下一步优先项**：① 修 `dispx` 只按全局轴的语义 ② 查上面那个反力方向异常
  ③ 换 `cart_patch.py` 的软骨网格（291 翻转 / 22 洞）

### 🧠 方法论教训（本轮第二次强化）

- **"文档说某配置不行"要复验**：`force` 的"病态"其实是**载荷被丢** + **关节张开**两件事叠加，
  两件都不是"刚体模态"。
- **能力已存在但要验语义**：`dispx`/`transverse` 早就在代码里，但一个按全局轴、一个
  从没和 `all` 做过对照 ⇒ **开关存在 ≠ 语义正确**。
- **读数出口是复现性的前提**：没有正式 reader，历史结论就无法复核（本条的实践依据）。
- **`rc` 要看 FEBio 的，不是 Python 的**：`ankle_lig_feb.py` 打印 `FEBio rc=`，
  而外层 `subprocess` 返回码恒为 0（本脚本已在 `acceptance.py` 里显式解析）。
""" .replace("LRN-NEW", f"LRN-{NEW}")

LT.write_text(t.rstrip() + ENTRY, encoding="utf-8")
print(f"LEARNINGS.md -> {LT.stat().st_size:,} B")

# TODO 更新
TD = ROOT / "ANKLE_GOAL_TODO.md"
s = TD.read_text(encoding="utf-8")
s = s.replace(
    "- [ ] **Step 1d** 写验收脚本 `scripts/ankle_fe/acceptance.py`；对当前配置跑出\"失败\" ⇒ 自证判据有区分度",
    f"- [x] **Step 1d** 验收脚本 `scripts/ankle_fe/acceptance.py` **已建成并跑通**（含 `.feb` 解析 + "
    f"Kabsch 转动拟合 + `<Loads>` 回归闸门）；`[LRN-{NEW}]`")
s = s.replace(
    "- [ ] **Step 2a** 韧带线软骨换 `cart_patch.py` 网格（`CART_DIR`）—— 先清掉 291 翻转/22 洞",
    "- [x] **Step 1e** 补**文本读数出口**（`LOG_DATA`，语法出 FEBio 手册 §3.19.1 + 附录 E.2.1）"
    f"⇒ `disp_*.txt`/`reac_*.txt`；`[LRN-{NEW}]`\n"
    "- [x] **Step 2-0** **找到稳定的姿势可用配置**：`DRIVER=dispx` + `FIX_TALUS=transverse` ⇒ "
    "NORMAL/5s/0 负 Jacobian，且与 `all` 响应不同（90 vs 162 N）⇒ **放开转动是可行的**\n"
    "- [ ] **Step 2a** 韧带线软骨换 `cart_patch.py` 网格（`CART_DIR`）—— 先清掉 291 翻转/22 洞\n"
    "  ⚠️ 格式不兼容：旧=`{j}_offsets`(hex8,8节点)，`cart_patch` 出=`{j}_nodes`+`{j}_tie_faces`(prism6,6节点)\n"
    "- [ ] **Step 2a2** 修 `dispx` 语义：当前只 prescribe `argmax(|LAX|)` 的**全局单轴**，非真轴向\n"
    "- [ ] **Step 2a3** 查清异常：`transverse` 反力在 x、`all` 反力在 z（顶面 BC 同为 prescribe x）")
TD.write_text(s, encoding="utf-8")
print(f"ANKLE_GOAL_TODO.md -> {TD.stat().st_size:,} B")

# session note 追加
SESS = WS / "memory/2026-10-04/ankle-fe-solve-and-handover.md"
SEC = f"""

---

## 十七、执行 `ANKLE_GOAL` Step 1~2（2026-10-05）→ `[LRN-{NEW}]` / `[ERR-20261005-007]`

### 17.1 三个真 bug / 真发现（按发现顺序）

1. **`<Loads>` 恒为空**（真 bug）：`ankle_lig_feb.py` 第 865~873 行，`continue` 之后整段载荷
   代码缩进在 `if` 体内 ⇒ **`DRIVER=force` 时载荷被静默丢弃**。修前 0 行 / 修后 6 行
   （|tr|=9.594 N/mm²）。`[ERR-20261005-007]`
2. **力控制失败的真因 = 关节初始张开**（不是"刚体模态"）：t=0 中位间隙 0.9732 mm、
   接触仅 31.0 mm²。`TOPFIX=lateral` 下胫骨轴向自由、接触几乎不承载 ⇒ 轴向近刚体模态
   + 病态接触。**证据**：`force` 在 **600 N 报 43 负 Jacobian（4 s）**、
   在 **0.6 N（1/1000）也不收敛（463 s）** ⇒ **与载荷大小无关**；而 `disp` NORMAL（15 s）。
   排除项：`FIX_TALUS=all` 无效、`FIB_SHARE=0` 无效、`NO_CONTACT=1` 仍报负 Jacobian。
3. **韧带线用旧软骨网格**：构建日志 `翻转修正 291 个；仍无法修正而剔除 22 个`、
   `det(J) min=5.08e-01`（`DET_MIN=0.5` ⇒ 勉强通过）。而 `build_feb.py` 线早已用
   `cart_patch.py` 修好。**且格式不兼容**：旧 `{{j}}_offsets`(hex8, 8 节点)
   vs `cart_patch` 出的 `{{j}}_nodes`+`{{j}}_tie_faces`(prism6, 6 节点)。

### 17.2 稳定的姿势可用配置（本轮最重要的正面结果）

| 配置 | 求解 | FEBio | 顶面反力 | \|u\|max | 负 J |
|---|---|---|---|---|---|
| `DRIVER=disp` | NORMAL ✓ | 15 s | — | — | 0 |
| **`dispx` + `transverse`** | **NORMAL ✓** | **5 s** | **[90.0,0,0] N** | **0.382 mm** | **0** |
| `dispx` + `all` | NORMAL ✓ | 31.7 s | [0,0,161.8] N | 0.535 mm | 0 |
| `force`（任意载荷） | ERROR ✗ | 4~463 s | — | — | 41~46 |

- **`transverse` 不是 no-op**（90 vs 162 N）⇒ 放松距骨 BC 真的改变行为
- **`dispx` 语义陷阱**：prescribe 的是 `argmax(|LAX|)` 指定的**全局单轴**（本例只 prescribe
  **全局 x**），**非真轴向** ⇒ 这解释 `[LRN-038]`「`dispx` 是 no-op」的印象
- ⚠️ **未解释异常**：`transverse` 反力在 x、`all` 反力在 z，而顶面 BC 同为 prescribe x。
  **在解释清楚前不要用这两个数字结算**

### 17.3 新增资产

- `scripts/ankle_fe/ankle_lig_feb.py`：加 `LOG_DATA` 文本读数出口
  （**语法出自 FEBio 手册 §3.19.1 + 附录 E.2.1，非猜测**）⇒ `disp_*.txt`/`reac_*.txt`
- `scripts/ankle_fe/acceptance.py`（7.1 KB，新增）：一条命令「跑 → 解析 → 出指标 → 判 C1」，
  含 `.feb` 解析、**Kabsch 刚体转动拟合**、`<Loads>` 空段回归闸门
- `temp/pyfebio_demo/run_lig.py`：Python 运行器（`env` dict 传参，绕开 cmd `set` 的坑）
- `temp/pyfebio_demo/scan_dead.py`：不可达代码扫描器（回测能抓住本 bug）

> **最关键的一条结构性教训**：脚本一直写着 `lig_read_{{tag}}.npz`「供 `read_lig.py` 用」，
> 但 **`read_lig.py` 从未存在**，`.xplt` 又是二进制 ⇒ **历史上所有读数都是临时脚本现算、
> 随 temp 清理丢失** ⇒ **这是"结论不可复现"的结构性原因**。现在有正式出口了。

### 17.4 与计划的差异（重要）

原 Step 2 = 「加弱弹簧 `TAL_SPRING`」。**实验后修正**：**换 `dispx` 已拿到稳定 + 放开转动，
不需要先实现弹簧**。前提是接受「位移控制 + 对 DISP 二分匹配目标力」（项目已有二分法传统）
⇒ **OpenSim 耦合仍可达**，且绕开整个力控制病态。
"""
if SESS.exists():
    tt = SESS.read_text(encoding="utf-8")
    if "十七、执行" not in tt:
        SESS.write_text(tt.rstrip() + SEC, encoding="utf-8")
    print(f"session note -> {SESS.stat().st_size:,} B")

D = WS / "memory/2026-10-04.md"
if D.exists():
    tt = D.read_text(encoding="utf-8")
    if "§十七" not in tt:
        D.write_text(tt.rstrip() + f"\n- `2026-10-04/ankle-fe-solve-and-handover.md` → §十七 执行 "
                      f"ANKLE_GOAL Step1~2（`[LRN-{NEW}]`/`[ERR-20261005-007]`）："
                      f"**`<Loads>` 恒为空是真 bug**（continue 后死代码）；力控制失败真因="
                      f"**关节初始张开**（0.6N 也失败 ⇒ 与载荷无关）；找到稳定姿势配置 "
                      f"**`dispx`+`transverse`**（NORMAL/5s，与 `all` 响应不同 90vs162N）；"
                      f"补文本读数出口 `LOG_DATA` + 验收脚本 `acceptance.py`\n", encoding="utf-8")
    print(f"daily note -> {D.stat().st_size:,} B")
print(f"⇒ LRN-{NEW}")
