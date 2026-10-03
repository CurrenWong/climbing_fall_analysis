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

