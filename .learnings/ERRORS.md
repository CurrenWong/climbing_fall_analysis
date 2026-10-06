# ERRORS — climbing_fall_analysis

> 命令/操作失败、静默失效、集成问题。改进建议见 `LEARNINGS.md`。

---

## [ERR-20261003-001] `t_max` 固定 2.0 s → h ≥ 19.6 m 静默返回 0 峰值

**类型**：静默失效（不报错、字段齐全、数值看似合理）
**严重度**：🔴 高 —— 越危险的场景输出越"安全"
**发现方式**：对标 [Y25] 扫 1–40 m 时看到 h≥21 m 峰值突然归零
**状态**：已修 + 已加守卫 + 已加回归测试

### 症状

```
h= 16.0  peak_f_pad=  122.20 kN
h= 18.0  peak_f_pad=  131.92 kN
h= 19.0  peak_f_pad=  136.63 kN
h= 20.0  peak_f_pad=    0.00 kN   ← 突然归零
h= 21.0  peak_f_pad=    0.00 kN
h= 25.0  peak_f_pad=    0.00 kN
```

**不抛异常，`max_compression_m` 也是 0，看起来像"这么高反而很安全"。**

### 根因

`climbing/pad.py` 的 `simulate_boulder_fall` 签名里 `t_max: float = 2.0`（固定）。

自由落体到接触需要 `t = √(2h/g)`：

| h | t_fall |
|---|---|
| 19.0 m | 1.968 s |
| **19.6 m** | **2.000 s** ← 分界 |
| 20.0 m | 2.019 s |
| 25.0 m | 2.258 s |

**h ≥ 19.6 m 时积分在撞击发生之前就结束**，全程 `comp = 0`，
于是 `in_contact.any()` 为 False → `i0 = i1 = 0` → 峰值取空切片 → 0。

### 为什么 v1 没发现

v1 的所有场景高度 ≤ 3 m（室内抱石墙 ~4.5 m），**永远碰不到这个边界**。
缺陷从 v1 起就存在，只是**参数空间没走到**。

### 复现

```python
from climbing.pad import simulate_boulder_fall
r = simulate_boulder_fall(height_m=25.0, mass_kg=77.0,
                          on_pad=False, posture="feet-first-stiff")
print(r.peak_force_n)   # 修复前: 0.0
```

显式给不够的 `t_max` 更能直接触发：

```python
simulate_boulder_fall(height_m=30.0, mass_kg=77.0, on_pad=False,
                      posture="feet-first-stiff", t_max=1.0)
# 修复后抛 RuntimeError: 积分结束时仍未发生接触 ...
```

### 修复

`src/climbing/pad.py`：

**① 自适应时长**
```python
t_fall = float(np.sqrt(2.0 * height_m / g))
if t_max is None:
    t_max = max(2.0, t_fall + 1.0)
```
签名从 `t_max: float = 2.0` 改为 `t_max: float | None = None`。
用 `max(2.0, ...)` 而非直接 `t_fall + 1.0`，是为了**不缩短**低高度场景的时长
（低高度下 1 m 的 t_fall+1 ≈ 1.45 s < 2.0 s，若直接替换会缩短，可能破坏既有行为）。

**② 守卫 —— 关键：让错误不可能静默通过**
```python
else:
    if height_m > 0.0:
        raise RuntimeError(
            f"积分结束时仍未发生接触：h={height_m:.1f} m 的自由落体耗时 "
            f"{np.sqrt(2.0 * height_m / g):.3f} s，而积分到 t={sol.t[-1]:.3f} s。"
            f"请增大 t_max（当前 {t_max:.3f} s）。"
        )
    i0, i1 = 0, 0
```

**只改默认值是不够的** —— 显式传入过小的 `t_max` 仍然会静默归零。
必须让"没发生撞击"这件事**抛错**。

### 回归测试（4 个）

`tests/test_climbing.py::TestIntegrationHorizon`
1. `test_high_fall_is_not_silently_zero` —— h=25 m 峰值必须非零
2. `test_peak_force_monotone_across_horizon` —— 跨 19.6 m 断点峰值仍单调上升
3. `test_insufficient_t_max_raises` —— 显式给 `t_max=1.0` 必须抛错
4. （在 `TestBoneSites` 里）高高度利用率 > 低高度利用率

### 更大的教训

见 `LEARNINGS.md` `[LRN-20261003-001]`：
**对标文献会要求你在更宽的参数范围上产生结果，从而逼出这类只在边界暴露的失效。**

---

## [ERR-20261003-002] 把 `f_flex` 同时喂给胫骨和腓骨（无载荷分配）

**类型**：物理建模错误（被回归测试抓出）
**严重度**：🟠 中 —— 只影响分部位排序，不影响总力
**状态**：已修

### 症状

新增的 `TestBoneSites.test_no_fracture_at_low_height_on_pad` 失败：

```
AssertionError: 屈膝缓冲（受控下跳） @ 1.0 m (软垫): 超过阈值 —— 腓骨两端(2.25)
  SiteLoad(key='fibula_ends', peak_force_n=12624.07, threshold_n=5600.0,
           utilization=2.254...)
```

**1 m 落在软垫上判腓骨骨折 2.25× —— 物理上荒谬。**

### 根因

`bone.py` 里 `peak_flex` 是**双腿合计**的轴向力（来自 `f_flex`），
我把它**原样**拿去比胫骨和腓骨各自的阈值。

但解剖上：**完整小腿的轴向载荷主要由胫骨承担，腓骨只分担约 10%**
（单侧 5%）。不分配 = 腓骨被高估约 20×。

### 修复

`BoneSite` 加 `load_fraction` 字段 + `load_from()` 方法：

```python
def load_from(self, peak_pad_n, peak_flex_n):
    gross = peak_pad_n if self.load_source == "pad" else peak_flex_n
    return gross * self.load_fraction
```

| 部位 | load_fraction | 理由 |
|---|---|---|
| 跟骨 | 0.50 | 双足均分 |
| 胫骨远端 / 骨干 | 0.45 | 双足均分 × 小腿内 ~90% |
| **腓骨两端** | **0.05** | **小腿内 ~10%，双足均分** |
| 股骨颈 | 0.50 | 双足均分 |
| 骨盆 | 1.00 | 单一结构 |
| 腰椎 | 1.00 | 单一柱 |

同一场景利用率 **2.25 → 0.11**。

### 教训

见 `LEARNINGS.md` `[LRN-20261003-003]`：
**凡"一个总力分配给多个并联构件"的地方，都要显式写下分配比例。**
比例写错是静默高估，且**只在某些参数点上暴露**（这里是 1 m 这种低位场景）。

**这个错误是我自己写测试时无意中抓到的** ——
再次印证 `阶段总结.md` §3.9："回归测试是最便宜的保险"。

---

## [ERR-20261003-003] 在诊断函数里硬编码了本该派生的数字

**类型**：自我错误（低成本，但同类反复出现）
**严重度**：🟡 低
**状态**：已修

### 症状

`bone.py::knockdown_report` 里我写死了结论文字：

```python
"结论：材料强度 × 承载截面 对跟骨**高估约 6–13 倍**。"
```

而同一个函数上面几行**运行时算出**的是：

```
所需结构性折减 = 16.9× … 33.8×
```

**同一段输出里，硬编码的 6–13× 和计算出的 16.9–33.8× 互相矛盾。**

### 修复

改成运行时派生：

```python
k_hi, k_lo = nominal / hi, nominal / lo
f"所需结构性折减 = {k_hi:.1f}× … {k_lo:.1f}×"
f"结论：材料强度 × 承载截面 对跟骨**高估约 {k_hi:.0f}–{k_lo:.0f} 倍**。"
```

### 教训

这正是 `阶段总结.md` §3.1 反复出现的那类错误
（**硬编码本该派生的值**）的又一次实例。
凡是"结论文字里带数字"的地方，都要问一句：**这个数字能不能算出来？**

---

## [ERR-20261003-004] 中文搜索工具 / 多行 `-c` 命令在 cmd.exe 下的老问题

**类型**：工具使用（非项目错误，但浪费时间）
**严重度**：🟡 低（每次浪费 1–2 分钟）
**状态**：规避

### 症状

用 `python -c "多行脚本"` 时反复报：

```
IndentationError: unexpected indent
```

原因是 cmd.exe 把换行/缩进吞掉或合并成一行。
另外 `python -X | tail -3` / `findstr` 在 cmd 下也常失败（无 `tail`）。

### 规避

- **多行 Python 一律写成临时 `.py` 文件再执行**，不要用 `-c`
- 需要看尾部输出时用 `more` 或直接重定向到文件再 `read_file`
- 用完即删（`_` 前缀标记为临时）

### 注意

这套环境**没有 pip**（`.venv\Scripts\python.exe -m pip` → `No module named pip`）。
装包用 **`D:\Program\uv\uv.exe pip install --python .venv\Scripts\python.exe <pkg>`**。

---

## [ERR-20261003-005] 解包字典时把「值」当成了「键」——`KeyError: '0.62'`

**类型**：实现错误（首次运行即暴露）
**严重度**：🟡 低（不静默，直接崩）
**状态**：已修 + 已加测试

### 症状

`ScenarioSampler` 初始化后第一次采样崩：

```
File "src/climbing/scenarios.py", line 489, in _mass
    _, mu, sd = MASS_MODEL[lab]
KeyError: '0.62'
```

### 根因

`MASS_MODEL` 是 `label -> (prob, mean, sd)`，我写的是：

```python
labels, p = zip(*[(v[0], v[1]) for v in MASS_MODEL.values()])
```

`zip(*pairs)` 把每对的**第一个元素**当 labels —— 也就是把 **prob 当成了
label、mean 当成了 prob**。于是 `rng.choice` 从 `[0.62, 0.36, 0.02]`
里抽，抽出浮点数 0.62 拿去查字典。

**正确写法**：

```python
self._mass_labels = list(MASS_MODEL.keys())
self._mass_p = np.array([v[0] for v in MASS_MODEL.values()])
```

### 教训

`zip(*[...])` 做转置式的解包**读起来不显眼**，很容易把元素顺序搞错。
当字典的 value 是元组时，直接 `list(d.keys())` + 列表推导取对应字段更清楚，
也更容易眼检。**别为了「一行」牺牲可读性** —— 这类错误不需要测试就能避免。

（对比：这个错误和 `[LRN-20261003-003]` 的载荷分配不同 ——
那个是**静默**高估、只在特定参数点暴露；这个是立刻崩，属于"好错误"。）

### 关联

同批还修了 `coverage_report()` 里两处**中文引号写成了 ASCII 双引号**
（`"踝"` 会被 Python 当成字符串边界，直接 SyntaxError）。
写中文文案时统一用 `「」`，别用 `""`。

---

## [ERR-20261003-006] pyfebio 0.3.0 能写不能读：`Literal[0,1]` 字段反序列化必炸

**类型**：第三方库 bug（往返失败）
**发现方式**：官方 README 样例跑完，用 `Model.from_xml()` 读回自己刚写的 .feb
**状态**：未修（上游 bug）—— 用法上按「只写不读」处理

### 症状

`Model.from_xml(feb_bytes)` 抛 77 个 ValidationError：
`boundary_.all_bcs.0.BCZeroDisplacement.x_dof: Input should be 0 or 1, input_type=str`

### 根因（最小复现见 `temp/pyfebio_demo/repro.py`）

pydantic 的 `Literal[0, 1]` **只接受 Python int，不接受 str `"1"`**：

| 输入 | `Literal[0,1]` | `int` |
|---|---|---|
| `1` (int) | OK | OK |
| `"1"` (str) | **FAIL** | OK |

而 pydantic-xml 反序列化时给的**永远是 str**。所以
`pyfebio/boundary.py` 里所有 `x_dof/y_dof/z_dof/relative/aggressiveness/
dtforce/cycle_buffer` 这些 `Literal[0,1]` 字段读回来必然失败。
`control_.solver.*.qn_method` 的 `type` 判别失败会连锁炸 20+ 个字段。

### 附带的调用坑

- `Model.from_xml(str)` → `ValueError: Unicode strings with encoding declaration
  are not supported`（.feb 带 `<?xml encoding='ISO-8859-1'?>`）→ **必须传 bytes**
- `Model.parse_file(path)` 也不行 —— 那是 pydantic 的 **JSON** 解析，直接
  `JSONDecodeError`（且已 deprecated）

### 正确读法

```python
m = pyfebio.model.Model.from_xml(Path("x.feb").read_bytes())  # bytes
```

但即使传了 bytes 仍会撞上面的 Literal bug。

### 结论

pyfebio 0.3.0 是 **generator-only** 的：写 .feb 可用，读回来不可用。
依赖它的流程只能是「Python 对象 → .feb → febio4」单向。
另外 `run_model()` 只是 `subprocess.run("febio4 ...")`，本机没装 FEBio
求解器时返回 rc=1 且只打印一行 cmd 错误，不抛异常。

**关联**：`temp/pyfebio_demo/{demo,roundtrip,repro}.py`

---

## [ERR-20261003-007] FEBio 4.13.0 不在系统 PATH；XPLT→HDF5 的 dataset 布局反直觉

**类型**：环境 + 集成坑
**状态**：已绕过（样例已跑通并通过物理校验）

### 坑 1：febio4 装在 D 盘但没进 PATH

安装位置 `D:\Program\FEBioStudio\bin\febio4.exe`（FEBio 4.13.0），
**不在系统 PATH 里**，也不在用户 PATH 里。每次跑都得手动：

```cmd
set "PATH=D:\Program\FEBioStudio\bin;%PATH%" && python demo.py --run
```

`run_model()` 是裸 `subprocess.run("febio4 -i ...")`，不做路径探测，
所以没 prepend 就静默返回 rc=1（**只打一行 cmd 错误，不抛异常**）。

### 坑 2：xplt.to_hdf5() 的 dataset 布局

`pyfebio.xplt` 需要**显式 import**（`__init__.py` 没导出），
转换入口是 **`xplt.to_hdf5(src, dst)`**（不是 `convert_xplt_to_hdf5`）。

转换后 layout（踩了两次才摸清）：

| 路径 | dtype | shape | 备注 |
|---|---|---|---|
| `meshes/0/nodes` | 结构化 `[id,x,y,z]` | (8,) | **是** recarray |
| `states/N/node_data/displacement/nodes` | `float32` | (8,3) | **不是**结构化，**列=xyz 分量** |
| `states/N/element_data/stress/box` | `float32` | (1,6) | Voigt 序 xx yy zz xy yz zx |

两个关键点：

1. **每个 state 是独立 dataset**（`states/0/...`、`states/1/...`），
   **没有**堆叠成时间轴的三维数组 → 不能 `f[key][:]` 一次拿全部帧。
2. 位移 dataset 是**纯 float32 (n,3)**，不是 recarray。
   写 `u[4]["z"]` 或 `u["z"][4]` 都会 `IndexError`
   （看着像结构化数组，习惯性写法必踩），正确写法 `u[4, 2]` / `u[:, 2]`。

### 顺带：`cmd.exe` 吞多行 Python

`python -c "..."` 里写多行 + 缩进，cmd.exe 会把换行折成空格导致
`IndentationError: unexpected indent`。**多行 Python 一律写成 `.py` 再跑。**

**关联**：`temp/pyfebio_demo/{demo,read_xplt,probe_h5}.py`

---

## [ERR-20261003-008] FEBio 4.x 手写 `.feb` 缺 `<solver>` 直接解析失败（且 `-silent` 掩盖原因）

**类型**：schema 漂移 + 诊断被掩盖
**严重度**：🟡 低（立刻失败，但报错难看到）
**状态**：已绕过（改用 pyfebio 生成）

### 症状

按记忆手写的"最小可解" hex8 模型，`febio4 -i min.feb`（verbose）报：

```
Reading file ...min.feb ...FAILED!
 *************************************************************************
 *                                ERROR                                  *
 * Component "" needs to have property "solver" defined (line 9)         *
 *************************************************************************
```

line 9 是 `<analysis>STATIC</analysis>` —— 明明没有叫 `solver` 的空组件，
实为 FEBio 4.x 要求 `<Control>` **必须**含 `<solver type="solid">` 子块。

同时踩到两个诊断坑：
1. 用 `-silent` 跑时控制台无报错、且**没有生成 `.log`** → 拿不到原因；
2. 真正权威的 schema 是 pyfebio 生成的 `my_model.feb`：`<Control>` 里
   `time_stepper` + `solver type="solid"` + `qn_method` + `linear_solver` 一整套。

### 根因

FEBio 4.x 的 Control 段比 2.x/3.x 复杂得多，手写极易漏字段；
且**解析阶段失败时 FEBio 不写 `.log`**，错误只出现在 stdout/stderr。

### 规避

- 工程/测试里**一律用 pyfebio 现建模型**（已验证 `IsotropicElastic` + 单 hex8 → `rc=0`），
  不手写 XML；
- 排查 FEBio 失败时**不要加 `-silent`**（会丢错误）；
- `climbing.coupling.febio_run.run_febio()` 默认 `silent=False`，失败时把
  stdout/stderr 尾部打进异常（并兜底找 `.log`）。

### 教训

"最小样例"如果**手写**，往往撞的是 **schema 坑**而不是物理坑。
能用官方库生成就别手写。

**关联**：`temp/pyfebio_demo/my_model.feb`（权威 schema）、
`src/climbing/coupling/febio_run.py`、`tests/test_febio_run.py`、`ERR-...-006/007`


---

## ERR-20261004-001 — FEBio `elastic-truss` 在真实混合模型里**读文件即崩**

**严重度**: 高（阻塞）
**现象**: 用 `<Elements type="line2">` + `<BeamDomain type="elastic-truss">`
（吴恺的"单轴梁单元"）加韧带，FEBio 4.13.0 **在 Reading file 阶段访问违例**
（returncode 3221225477），无任何 ERROR 文本。

**排查（做了完整隔离）**:
- 最小实验 **全部通过**：hex8+杆、tet4+杆、杆两端共享实体节点、杆一端独立、
  2 根杆、+`<Surface>`+`<SurfacePair>`+`<Contact>` — 都 `NORMAL TERMINATION`。
- 真实模型（32k 实体）：**逐项减法**定位到
  「删掉 `<Elements type="line2">` 块 ⇒ 立刻能读」（`04_no_truss` READ-OK）。
- 换单条韧带（`LIG_SUB=1..10`）**每条都崩** ⇒ 不是某条韧带。
- 换 `type`：`linear-truss` 拒绝标准材料（`invalid value for attribute "type"`）；
  `linear-beam` 也崩；`<area>` 不是合法材料标签。

**结论**: `elastic-truss` 是手册标注的 **EXPERIMENTAL** 特性，在真实尺度混合
网格上不可靠。**不要在它上面继续投入**。

**替代路线**（按把握排序，见 `ANKLE_LIG_TODO.md` 第二节）：
① 把韧带专属节点焊到最近骨表面节点，用 shell；
② 重做骨网格保留全部 THUMS 表面节点，使韧带与骨天然共节点。

---

## ERR-20261004-002 — FEBio shell 与 solid 共用节点：**只属壳的节点转动自由度无约束**

**严重度**: 高
**现象**: 用 THUMS 韧带的原始 quad4 网格做成 `ShellDomain`，文件能读、能进求解，
但**首个载荷步就报 22 个 negative jacobian**，在 **1/220 的载荷**（time 0.00454545）
就失败。

**根因**: 韧带引入了 **137 个「只属韧带」的节点**（11,803 − 11,666）。
这些节点只被壳单元用到、不被任何实体用到 ⇒ 壳的**转动自由度（含钻孔转动）
没有约束** ⇒ 刚体模态 ⇒ 刚度阵奇异 ⇒ 位移解无意义 ⇒ 单元翻负。
**判据**：失败发生在 **1/220 载荷**、且**位移量级完全不合理**时，
优先怀疑**奇异/刚体模态**，而不是"材料太软压溃了"。

**尝试过的修法及结果**:
- 给这些节点加转动零位移 BC：`<u_rotation>` **不是合法标签**
  （FEBio 手册里 `*_rotation` 只有**刚体**转动约束）；本轮未找到壳节点转动的 BC 写法。
- 改成"只用骨节点"的 quad 条带 ⇒ 节点数确实回到 11,666（**零新增** ✓），
  但 truss 配对共享节点 ⇒ 只生成 **12** 个合法 quad（太少）。
- 改成"只用骨节点"的 tri3 三角网（42 个）⇒ **读文件即崩**（tri3 作壳域不支持）。

**下一步**: 焊节点（方案①）或重做骨网格（方案②）。

---

## ERR-20261004-003 — 空 `<Elements>` 块 / `Shelldomain type="default"`

两处 FEBio 的"报错信息与病因不匹配"，各花掉一轮：
1. **元素数为 0 的 `<Elements type="X" name="Y"></Elements>`**
   ⇒ 报 `tag "Elements" (line N): missing attribute "id"`（完全没提"空块"）。
   修法：空块**不写**。
2. **`<ShellDomain ... type="default">`** ⇒ 报
   `invalid value for attribute "default"`（把值当成了属性名）。
   修法：`type` 是**可选**的，**省略**即可（FEBio 按材料/单元类型自行推断）。


---

## ERR-20261004-004 — 求解器配置问题被误当成物理问题（"每步 97 次迭代磨不完"）

**严重度**: 中（浪费大量时间）
**现象**: 加了韧带的踝关节模型在力控制下：
- 600 N 时 **首个载荷步（1/220 载荷）** 就报 22 个 negative jacobian；
- 降到 6 N 能跑 3 步，但报 **"Max nr of iterations / Max nr of reformations reached"**；
- 每步 **~97 次牛顿迭代**，3 步花 **407 s**。

**误判过程**（我的错误链条）：
1. 以为是**刚体模态**（LRN-016/017 的教训）⇒ 去固定舟骨、去写转动约束 BC —— 都没用
2. 以为**韧带壳与骨共用节点**引入自由转动自由度 ⇒ 改成只用骨节点的条带/三角网 —— 也没用
3. 以为 **`elastic-truss`** 是正解 ⇒ 真实模型读文件就崩（见 ERR-20261004-001）

**真因**（读 FEBio 迭代日志才看出来）：
```
convergence norms :     INITIAL         CURRENT         REQUIRED
   residual            9.840343e+06    3.239042e+01    0.000000e+00   ← 要求降为 0，永不可能
   displacement        1.870091e-01    4.064515e-02    1.257888e-07   ← 容差过严，差 5 个数量级
```
- **残差判据 REQUIRED = 0** ⇒ 该判据形同废除（**同一个坑 MEMORY 里已记过一次**）
- 实际只能靠位移判据，而它 **1.26e-07 过严** ⇒ 死磨
- 缺 `<qn_method type="BFGS">`（FEBio 的标准拟牛顿加速）

**修法**：按手册 §3.3 写全 `<solver>`：`max_refs/diverge_reform/reform_each_time_step/
dtol/etol/rtol/lstol/min_residual/rhoi/qn_method(BFGS)` + `<time_stepper>` 的
`dtmin/dtmax/max_retries/opt_iter`。

**教训**：
1. **"失败在 1/220 载荷"+"每步上百次迭代" ⇒ 先查求解器配置，别先改物理/几何**
2. **必须读 FEBio 的迭代日志**（`convergence norms` 那几行），它直接告诉你哪条判据没满足；
   只看"ERROR TERMINATION"会盲猜
3. **`REQUIRED = 0` 是个红旗**：任何"要求降到 0"的收敛判据都不可能满足

---

## ERR-20261004-005 — pyfebio 读 shell/beam 域的 .xplt：**两处** KeyError

**严重度**: 低（有补丁）
**现象**: `pyfebio.xplt.to_hdf5()` 对含 shell/beam 域的模型失败：
1. `parse_mesh` → `KeyError: 'name'`
   （`domain["name"]` —— 注意是**双引号形式**，只替换 f-string 里的 `{domain['name']}`
   **不会命中**，本轮因此白试一轮）
2. `parse_state` → `KeyError: np.int32(4)`
   （`mesh_dict[set_lut[key]][set_id]` —— **`set_id` 与 `set_id + 1` 两处都要兜底**）

**修法**（`temp/pyfebio_demo/read_lig.py` 已实现）：用 `inspect.getsource` 取源码、
字符串替换加 `.get(...)` 兜底、`exec` 回写模块。

**附带**：`<Nodes>` 的坐标是**结构化数组** `(id, x, y, z)`，不是 `(n,3)`；
要 `np.column_stack([c["x"], c["y"], c["z"]])`。

---

## [ERR-20261005-007] `ankle_lig_feb.py` 的 `<Loads>` 段**恒为空** —— `continue` 之后的载荷代码全是死代码

- **日期**：2026-10-05
- **文件**：`scripts/ankle_fe/ankle_lig_feb.py` 第 865~873 行
- **严重度**：🔴 高（静默、且导致错误结论进入文档）

### 症状
`DRIVER=force`（力控制）跑不动：600 N 时报 43~46 个 negative Jacobians，首步即失败。
文档把它记为「**力控制病态（胫骨只靠接触约束 ⇒ 刚体模态）**」，并据此改用 `DRIVER=disp`。

### 真因
```python
for surf, F, _faces in (("tibia_top", F_t, top_tibia), ...):
    if not _use_force or not _faces or F <= 0:
        continue
        # traction = 力 / 面积 ...        ← 缩进 8，在 if 体内
        a = area_proj(...)                 ← 以下 7 行全是死代码
        w(f'    <surface_load ...')        ← 永不执行
```
**`DRIVER=force` 时整段 `<surface_load>` 从未写出 ⇒ 模型零载荷。**
（无载荷 + 顶面轴向自由 ⇒ FEBio 空磨/奇异。**"刚体模态"只是零载荷的后果，不是根因。**）

### 实测证据
| | `<Loads>` 段 |
|---|---|
| 修前 | **0 行**（空） |
| 修后 | 6 行，`traction=-7.821332,-1.294633,-5.402078` ×2 ⇒ \|tr\|=**9.594 N/mm²** |

旁证：`temp/pyfebio_demo/mix_probe/B_one_free.log` 里 FEBio 自己的警告
`* No force acting on the system. *`

### 修法
把第 867~873 行 dedent 4 空格到 `for` 体内（`edit_file` 已应用；`ast.parse` OK）。
已加注释标注本 ERR 编号。

### 复发防线
1. `temp/pyfebio_demo/scan_dead.py` —— **不可达代码扫描器**（AST 判 `body[i+1]` 是否在
   `return/continue/break/raise` 之后）。**回测确认能抓住本模式**（对修前片段报 1 处）。
   实测 `ankle_lig_feb.py` / `build_feb.py` / `cart_patch.py` / `cart_refine.py` / `acceptance.py`
   现均 **0 处**。
2. `scripts/ankle_fe/acceptance.py` 内置 **`<Loads>` 行数回归闸门**（`DRIVER=force` 时行数须 > 0）。

### 教训（已并入 `[LRN-20261005-058]` 的方法论）
- **凡"某驱动方式不行"的结论，先验产物（`.feb` 的载荷段/BC 段）再看求解器报错。**
  光看求解器报错永远看不出"载荷根本没施加"。
- **缩进型静默 bug 是 Python 写代码生成器的高危区**：`continue`/`return` 之后多缩进一层，
  逻辑不报错、语法不报错、结果全错。

---

## [ERR-20261005-008] FEBio `<traction lc="...">` 被静默忽略 —— 载荷曲线必须挂在 `<scale>`

- **症状**：多步分析 step2（力控制）**无论总载荷多大**（100 / 600 / 5000 N 实测），
  都在第一个增量报 42~48 个负 Jacobian，且**看似**"卡在 11~12 N"。
  更早的表现：`DRIVER=force` 历史上一律失败，看起来像"力控制病态"。
- **根因**：写成 `<traction lc="2">-0,-0,-9.93</traction>`。
  **FEBio 不认 `<traction>` 上的 `lc`，静默忽略 ⇒ 满载荷在第一增量瞬加。**
- **正解**（手册 §3.13.2.2 / p121 原文示例）：
  ```xml
  <surface_load type="traction" surface="tibia_top">
    <scale lc="2">1.0</scale>      <!-- ★ 载荷曲线挂在 scale 上 -->
    <traction>0,0,-9.93</traction>
  </surface_load>
  ```
  原文："An optional load curve can be defined for the **scale** element using the `lc` attribute."
- **诊断信号**（有用）：**失败点对总载荷量级不敏感** ⇒ 一定是"满值瞬加"而不是"斜坡"。
  用 `载荷 = FORCE_N × (t-1)` 反推会得到假象的"绝对值天花板"。
- **影响面**：此前所有 `DRIVER=force` 实验都未真正测过"斜坡力控制"。
- **修法**：`ankle_lig_feb.py` 的 traction 输出改为 `<scale lc=..>1.0</scale>` + `<traction>`；
  单步路径同步加 `<scale lc="1">1.0</scale>`。

---

## [ERR-20261005-009] 加载方向符号反了 —— 整条踝线一直把关节**拉开**（拉伸）而非压实

- **症状**（长期误判为"力控制病态"）：`DRIVER=force` 无论如何调参都在 ~10-12 N 崩掉；
  单元翻负；换弹簧刚度/接触算法/步长/落位量全无效；结果不可复现（同一配置 337 N vs 159 N）。
- **根因**：载荷与落位位移都用了 `-LAX`（LAX = 胫骨→距骨 远端方向；ALIGN_AXIS=1 后 = +z）。
  **实测证明 `-LAX` 是"把胫骨往近端拉"= 把关节拉开**，抵抗者是韧带
  ⇒ 能承受的量级正好是 ~10-200 N ⇒ 一切"卡住"都发生在拉伸侧。
- **决定性证据**（都收敛，符号相反）：
  | DISP_MM | NORMAL | 顶面 Fz |
  |---|---|---|
  | **−0.25** | ✅ | **−132.56 N（压缩）** |
  | **+0.25** | ✅ | **+187.77 N（拉伸）** |
  历史上记录的 `187.77 / 187.90 / 202.29 / 209.5589 / 592.61 N` **全部是 +DISP 的拉伸值**。
- **修正**：新增 `LOAD_SIGN`（默认 -1 保持旧行为），压缩方向 = `+LAX`（`LOAD_SIGN=+1`）；
  `_topdisp()` 的落位位移同步改为压缩方向。
- **修正后结果**（`dispx`、压缩、PENALTY、全程 NORMAL）：
  | DISP_MM | Fz (N) | k (N/mm) |
  |---|---|---|
  | −0.25 | −132.56 | 530 |
  | −0.50 | −292.22 | 584 |
  | −0.80 | −479.60 | 600 |
  | **−1.00** | **−608.82** | **609** |
  | −1.20 | −741.55 | 618 |
  | −1.60 | −1017.36 | 636 |
  ⇒ **600 N 压缩态已被干净达到**（且余量到 1017 N）。
- **教训**：**"力控制病态"的排查应先验证受力方向的物理正确性**（压缩 vs 拉伸），
  再怀疑算法。一个符号可以让所有参数调优都变成噪声。
- **连带更正**：`[LRN-064]` 的"1.60 mm → 592.61 N（窗口 2.5×）"是**拉伸**结果，
  不能再当作"600 N 到不了被推翻"的证据；压缩侧本来就通。

---

## [ERR-20261005-010] `_outer()` 只认旧软骨键 ⇒ `POSE_TILT` 一直崩；且被"旧 `.feb` 陷阱"掩盖

- **症状 A**：任何带 `POSE_TILT` / `POSE_TILT_ML` 的运行 `rc=1`，
  `KeyError: 'tibiotalar_tibia_offsets is not a file in the archive'`。
- **根因 A**：`_outer()`（姿势的关节面间距自检用）只读**旧格式**键 `{j}_offsets`；
  而 `CART_NEW`（`cart_patch.py`）用 `{j}_nodes` + `{j}_is_offset`（节点键名 `cn:{j}:{k}`）。
  ⇒ **姿势功能自新软骨格式起就是坏的**（`[LRN-20261005-063]` 之后的所有姿势结论都不可靠）。
  修法：`CART_NEW` 分支走 `cn:{j}:{k}` + `CN_ISOFF[j][k]`。
- **症状 B（更隐蔽）**：8 个不同 `POSE_TILT` 跑出**完全相同**的 `Fz = −608.82 N`、剪切恒 0。
- **根因 B**：脚本没查**构建返回码**、也没**删旧 `.feb`**。构建失败后
  **FEBio 会读上一次留下的 `.feb`** ⇒ 所有姿势其实跑的是同一个文件。
- **教训（两条都要写进 runner 模板）**：
  1. **构建完必须查 `rc`，并确认 `.feb` 的 mtime 是本次的**；
  2. **每次跑前删旧 `.feb` 与旧读数** —— "结果逐字相同"永远是第一嫌疑。
- **另一条同源教训**：**反力必须看符号**（压缩为负）。只看幅值会把**拉伸**当成成功
  （`[ERR-20261005-009]`）。

---

## [ERR-20261005-011] 几何旋转式姿势的"落位"符号反了 → 姿势功能被误判"不可用"

- **症状**：`POSE_DF` / `POSE_IE`（绕踝中心刚性旋转足部组的**真实内外翻/屈伸**姿势）
  一律在 t=0 报 `Negative jacobian`，`DISP_MM` 从 0.25 到 1.2 全 FAIL
  ⇒ 此前结论"几何旋转式受阻"（写进了使用说明 §3）。
- **根因**：旋转后要沿加载轴平移足部组把关节面"落位"回去，原实现是
  `_shift = _d1 - _d0`（`_d1` 转后间距、`_d0` 转前），然后沿 +axis 平移。
  **实测 `POSE_IE=+10°`：转前 0.0211 → 转后 0.3625 → 平移 +0.3414 后 0.4416**
  ⇒ **平移把间距越推越大（符号反）** ⇒ 关节在 t=0 就张开 0.44 mm ⇒ 负 Jacobian。
- **修复过程中的两个错误尝试（都记下来，避免重犯）**：
  1. **割线法**（2 次测量解 `gap(s)=d0`）⇒ 残差仍 +0.33 mm。原因：**"最近点对"随平移换点**，
     `gap(s)` 非光滑，一步割线失准。
  2. **括入 + 二分**（假设 `gap(s)` 单调）⇒ 更差（0.4664）。原因：
     **`gap(s)` 是 V 形（有唯一极小），不是单调** ⇒ 对目标值二分无解。
  3. **无界一维极小化** ⇒ 给出 +2.04 mm 的荒谬平移（关节面被滑到无关位形、韧带预拉伸）
     ⇒ **平移量必须设界**。
- **正解**：对平移量 s 做**有界（|s| ≤ `POSE_SEAT_LIM`=0.6 mm）一维极小化**
  （粗扫 0.01 mm + 细化 0.001 mm），目标改为**可达的最小间距**（Anderson 的 minimal contact），
  而不是"回到转前的某个值"。
- **修后结果**：`POSE_IE` = ±1/±2/±3/±5° **全部 NORMAL**（落位后间距 0.03~0.15 mm）；
  `POSE_DF=2°` NORMAL、`5°` 仍失败 ⇒ **可用范围 ≈ ±5°**。
- **教训**：① "落位/贴合的平移方向"必须**实测验证**（看间距是变大还是变小），不能靠推导；
  ② 非光滑目标函数不要用割线/二分；③ **优化变量要设物理界**，否则会被"假最优"带走。

## ERR-20261006-001 — `extreme_faces()` 选出的 NodeSet **跨姿势不稳定** ⇒ 跨姿势比较会错

**症状**：用同一个 `NodeSet` 名从中性网格与姿势网格各取坐标做差，报
`operands could not be broadcast together with shapes (227,3) (194,3)`。

**根因**：`tibia_top_ns` / `calcaneus_bottom_ns` / `talus_low_ns` 由 `extreme_faces()`
按「带内极值面 + 法向判据」选取 ⇒ **姿势一变，符合条件的节点集合就变**
（实测 `calcaneus_bottom_ns`：中性 **194** → 某姿势 **227** 节点）。

**为什么危险（比报错更坏的情形）**：集合大小偶尔相同的情况下不会报错，
**却是在比较不同物质点** ⇒ 结论静默错。

**正解**：跨姿势/跨运行比较**一律固定用其中一个构型（通常 r=0）的 ID 列表**，
按同一批 ID 从两边取坐标。`check_c2.py` 已按此实现（`NS_CAL/NS_TAL/NS_TIB`）。

---

## ERR-20261006-002 — 我的"`POSE_IE` 标签反了"结论是错的；根因是**循环论证**

**错的结论**（曾写入 `MEMORY.md` / 使用说明 / `ANKLE_FIX_TODO.md`）：
「`POSE_IE` 打印"内翻+"是反的，**+IE 实际是外翻**」。

**它是怎么产生的（循环论证链）**：
1. 观察到「CFL 应变在 `POSE_IE=−3` 时最高」
2. 由此**反推**「所以 −3 必是内翻（因为内翻才拉外侧韧带）」← **循环：用数据本身定标签**
3. 再把该标签回代，宣布「内翻拉外侧 ⇒ 拮抗关系物理正确」← **用自己推出的前提验证自己**
（`temp/pyfebio_demo/save_atfl.py` 里 "现象：内翻（`POSE_IE=−3`）…" 是这链条的化石）

**打破循环的证据**（`check_c2.py`，两个**不依赖模型符号约定**的解剖锚定判据）：
- `POSE_IE=+5` ⇒ 足底法向**外侧分量减小** −0.0825 且**外侧缘下沉** −0.81 mm（内侧 +0.82）
- ⇒ **+IE = 内翻**（与代码注释一致）；CFL/ATFL/PTFL **长度 +2.47/+1.75/+0.63 mm** 同步佐证

**真正的错在哪**：不是标签、也不是应变算错，而是**引用了错的口径** ——
旧结论引用的是**(b) 载荷增量**，却当成总应变量解读。详见 `[LRN-20261006-075]`。

**教训**：① 凡"用数据反推标签、再用标签证明数据"的链条，必须找一个**外部锚定**
（这里是斐骨=外侧、足底法向）来打破；② 同一物理量有多个口径时，**先写清口径再比较**。

## ERR-20261006-003 — **我上轮的"加载方向符号修正"把方向搞反了**（`ERR-20261005-009` 本身是错的）

**现状**：`DRIVER=dispx` 下 **`DISP_MM > 0` 才是压实（推向关节）**，
`DISP_MM < 0` 是把关节**拉开**。2026-10-05 我判定"正解是 `LOAD_SIGN=+1` + `DISP_MM<0`"，
**方向判反了**，并把本来正确的历史约定（`DISP_MM>0`）标成"拉伸、作废"。

**三条**互相独立**的证据（都可复现）**
1. **`.feb` 实体**：`DISP_MM=+1.0` ⇒ 写入 `<value>−1.000000</value>`；`−1.0` ⇒ `+1.000000`
   ⇒ **z 值 = −DISP_MM**
2. **几何（不依赖任何符号约定）**：顶面 z̄=−653.5、距骨 z̄=−709.4（**距骨在下方**）；
   `DISP=+1.0` ⇒ 顶面 dz=**−0.927**（向下、朝距骨）⇒ **靠近**；
   `DISP=−1.0` ⇒ dz=**+1.000**（向上、离开）⇒ **远离**
3. **接触力（距骨全固定、无韧带）**：`DISP=+1.0` ⇒ **149.59 N**（接触咬合）；
   `DISP=−1.0` ⇒ **0.0000 N**（接触分离）
   —— 反向时若有韧带，力全部来自**韧带拉伸**（`LIG_MODE=none` ⇒ 0 N）

**旁证（闭环）**：历史零回归值 `TIE=1 CART_NEW=1 DISP_MM=+1.6` ⇒ **592.6020 N** 逐位复现
⇒ **历史约定本来就是压实且可复现**；我把它标"拉伸/作废"是错的。
`LOAD_SIGN` 在 `dispx` 分支**不改变写入的 z 值**（实测 4 组：z 只由 DISP 的符号决定）。

**修正**：压痛方向一律用 **`DISP_MM > 0`**；`LOAD_SIGN` 对 `dispx` 无效，
不要再靠它"纠符号"。**本文件里 `ERR-20261005-009` 的"正解"部分已作废。**

---

## ERR-20261006-004 — `ligaments.py` 两处真 bug（但**只影响 `LIG_MODE=truss`**）

1. **纤维轴取错**：原用「整片点云主成分最大奇异向量」，对**宽而短**的韧带与骨-骨方向差
   **PTCL 89.3°、SYND 84.0°** ⇒ `L` 从 18.0/32.4 被抬到 50.7/44.5 ⇒ `A_cross=面积/L×厚` 偏小
   ⇒ **杆模式下轴向刚度偏低 PTCL 2.81× / SYND 1.37×**（长条型不受影响：CFL 0.7°、ATFL 2.1°、
   MD1~4 1.3~4.8° ⇒ 0.98~1.09×）。已改为**骨-骨方向**。
2. **配对退化成扇形**：原「按沿轴归一化位置最近」在 ATFL 上让 **4 根杆全部汇聚到同一腓骨节点**
   ⇒ 足迹被压成一点。已改为**1:1 贪心 + 剩余端点各自就近**（实测 10 条韧带足迹覆盖 **100%**）。
3. **⚠️ 影响范围**：`ankle_lig_feb.py` 里 `a_cross` 只用于 `LIG_MODE=truss`；
   默认 **`LIG_MODE=shell` 用固定 `E=260`** ⇒ **默认路径的 FE 结果不受影响**（实测
   `DISP=+1.6` 修复前后均 592.6020，逐位一致）。应变分析走杆对，故数字有微变（CFL 11.79→11.48）。


---

## [ERR-20261006-005] 加载面积 `top_area` 与几何不一致（11~45×）⇒ 力控模式结果不可信

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
