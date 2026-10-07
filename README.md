# 室内抱石与攀岩坠落损伤的多层级数值建模

### climbing_fall_analysis

**Multi-level numerical modelling of indoor bouldering and climbing fall injury.**

室内抱石坠落是低能量、高发生率的损伤场景，急诊中半数以上伤及下肢，踝损伤居首。这个仓库回答一个问题：软垫到底保住了什么、又牺牲了什么。它由一条 1D 软垫-人体冲击模型起步，逐步扩展出 2D 转动、踝骨与韧带双判据、头颈惯性分支，以及一条把地面反力（GRF）串到 OpenSim 多体动力学和单骨有限元的载荷链，最后用蒙特卡洛把临床流行病学权重带进模型输出。数值按文献靶子做验证，结论按可信边界分层表述。许可为 MIT。

> 项目的第一原则是诚实：**相对排序与机制可引用，绝对量值与骨折率不可引用**。请在引用任何数字前先读[已知局限](#已知局限-important)。

---

## 核心结论 Key Findings

- **软垫护「压」、伤「扭」。** 同高 2 m 对照：软垫峰值 34.8 kN / 接触 34 ms，刚性地面 42.4 kN / 20 ms（峰值力约 -18%，接触时长约 +70%）。削峰、摊时间的收益是轴向的；代价是足部沉入增大踝旋后角，旋转失效风险上升。
- **踝损伤首次进入模型输出。** 1000 次蒙特卡洛给出踝损伤 34.5%（对临床实测 40.0%，1.16×），踝扭伤 34.5%（对 28.0%，1.23×）。此前 1D 纯压缩模型结构上算不出扭伤。
- **复现了文献的两条反直觉结论。** 同角度下内翻损伤轻于外翻（力臂随角度符号切换）；足在垫上的滑动与卡住是几何相变，只由落地角与滑动摩擦系数的相对大小决定，与坠落高度无关。
- **绝对量不可移植。** 独立 FE 线与 1D 线的首折高度差 3-4 倍；刚性地面文献（冲击 1-5 ms）与软垫工况（冲击窗 56-67 ms）分属两种情形，绝对峰值时刻与骨折率不可直接迁移。
- **验证靶子分层。** V2 / V4 / V6 / V9 达成，V8 的单调性达成；V1 部分达成；V3 / V5 / V7 / V8 的绝对值未达成。根因归为两类：工况错配，或模型范围局限。

---

## 安装与快速开始 Quick Start

### 环境要求

来自 `pyproject.toml`：

- `requires-python = ">=3.11"`
- 运行依赖：`numpy>=2.0`、`scipy>=1.13`、`matplotlib>=3.8`、`pandas>=2.2`
- 测试依赖（可选）：`pytest>=8.0`（`[project.optional-dependencies] dev`）
- 构建后端：`hatchling`；发行名 `climbing-fall-analysis`，导入名 `climbing`

### 安装

```bash
pip install -e .            # 安装本项目
pip install -e ".[dev]"     # 连同 pytest 一起装（跑回归用）
```

若不想安装，用 `PYTHONPATH` 直接指向 `src`（PowerShell）：

```powershell
$env:PYTHONPATH="src"
python -c "from climbing.pad import simulate_boulder_fall; print(simulate_boulder_fall().peak_force_kn)"
```

`pyproject.toml` 已写入 `pythonpath = ["src"]`，装了 pytest 后也能直接找到包。

### 最小可运行示例

```python
from climbing.pad import simulate_boulder_fall

# 2 m 自由落体，77 kg，受控屈膝落地，默认 20 cm PU 软垫
res = simulate_boulder_fall(height_m=2.0, mass_kg=77.0, posture="controlled-drop")

print(f"峰值垫反力: {res.peak_force_kn:.1f} kN")
print(f"最大压缩:   {res.max_compression_m * 100:.1f} cm")
print(f"是否压穿:   {res.bottomed_out}")
print(f"残余能量:   {res.energy_residual_j:.1f} J")
```

`POSTURES` 内置 10 种落地姿势（如 `flat-flop`、`butt-impact`、`head-first`、`toe-point`、`tuck-roll`）。2D 版本接口同形：

```python
from climbing.pad2d import simulate_boulder_fall_2d

res2d = simulate_boulder_fall_2d(height_m=2.0, mass_kg=77.0)
```

### 复现主要结果

| 命令 | 作用 |
| --- | --- |
| `python scripts/phase_p1_montecarlo.py` | 场景库蒙特卡洛，n=1000，seed 20261003，输出 `results/p1_montecarlo.txt` |
| `python scripts/validate_vs_fem.py` | 对标文献的分部位阈值验证 |
| `python -m pytest tests -q` | 全套回归（626 passed / 8 skipped / 0 failed） |

脚本输出写到 `results/`，图表写到 `results/figures/`。

`scripts/opensim_fe/` 与 `scripts/ankle_fe/` 是载荷链与踝 FE 的复现脚本，依赖 `model/` 下的第三方模型与 FEBio / OpenSim 运行时，**不装这些无法运行**（见[第三方模型与文献](#第三方模型与文献重要)）。

---

## 目录结构 Repository Layout

```
climbing_fall_analysis/
├── README.md              # 本文件
├── LICENSE                # MIT
├── CITATION.cff           # 引用元数据（作者字段为占位）
├── pyproject.toml         # 打包与依赖（Python >= 3.11，hatchling）
├── src/climbing/
│   ├── pad.py             # 1D 软垫-人体冲击：两质量单边接触 + 泡沫本构 + 腿屈曲
│   ├── pad2d.py           # 2D 扩展（z, x, θ）+ 踝姿角 + 滑动/卡住相变
│   ├── ankle_injury.py    # 踝骨利用率 + 韧带扭伤 / 撕脱骨折判据
│   ├── head_branch.py     # 头-颈惯性分支（定高坠落峰值时刻）
│   ├── bone.py            # 7 个骨部位的失效阈值
│   ├── scenarios.py       # 数据驱动场景库 + 覆盖率报告
│   ├── injury.py          # 分部位损伤评估
│   ├── metrics.py         # HIC、峰值 g 分档、风险汇总
│   ├── rope.py, belay.py  # 绳索坠落通道（Phase A，保留）
│   ├── bodies.py          # 多刚体链（结构完成，未通过验收，暂不用）
│   └── coupling/          # GRF -> OpenSim -> 单骨 FE 的载荷链
├── scripts/
│   ├── phase_p1_montecarlo.py   # 主结果：场景库蒙特卡洛（n=1000）
│   ├── opensim_fe/              # 多体 + 单骨 FE 复现脚本
│   └── ankle_fe/                # 踝 FE 建模脚本（THUMS / 韧带）
├── tests/                 # 22 个测试文件（回归套件）
├── results/               # figures/ + opensim_fe/，全部可再生
├── docs/                  # 论文稿、方法、报告、S5 接触系列、工作台账
│   ├── README.md          # 文档索引（建议从这里进入）
│   ├── paper/             # 论文初稿 / v2 数学建模版 / 出处索引 / 评审 / 汇总
│   ├── methods/           # 方案更新 v2、各扩展方案与建模交接
│   ├── reports/           # 阶段结论、进度、对齐与文献解析（含 S5_QUOTABLE）
│   ├── s5-contact/        # S5 脚-垫接触收敛系列（17 篇）
│   └── notes/             # 阶段总结 + 24 个工作台账
└── legacy/                # 早期 1D 代码，仅供追溯，勿引用其数值
```

**不随仓库分发**：`model/`（约 320 MB 第三方模型）、`paper/`（约 163 MB 文献 PDF）、`temp/`（约 6 GB 仿真产物，可再生）。三者均在 `.gitignore` 中，原因见[第三方模型与文献](#第三方模型与文献重要)。

---

## 模型概览 Model Overview

各层的数学表述与推导见 [`docs/paper/论文v2_数学建模版.md`](docs/paper/论文v2_数学建模版.md)，本节不含公式。

- **1D 软垫-人体冲击。** 把人体压成两个质点，下肢与躯干用屈曲弹簧串联，底端与软垫发生单边接触。软垫用开孔 PU 泡沫的「平台 + 压实」应力-应变本构，接触面积随沉入渐进展开。人体关节屈曲是第二个吸能器，它与软垫的串联行程共同决定峰值力。
- **2D 扩展。** 在垂直位移之外加入水平位移与转角，让落地姿角与足部沉入可被表达。踝姿角驱动旋转派生量，并带来滑动与卡住的几何相变，内翻与外翻通过力臂符号切换。
- **失效判据。** 骨侧用弯曲加横向剪切算利用率，韧带侧判扭伤，另有撕脱骨折判据，三者共同给出踝的结论。
- **头-颈惯性分支。** 以躯干轨迹为基础激励，把头颈当作受迫振子，给出头部应力峰值时刻，用于对标文献的时序靶子。
- **载荷链。** 地面反力经 OpenSim 多体动力学得到关节反力，再作为单骨有限元（FEBio）的边界条件，得到应力场。这是把宏观运动学接到组织级应力的通道。
- **蒙特卡洛。** 在有效包线内按临床频次抽样 1000 次（seed 20261003），把姿势、高度、质量等自变量组合成场景，汇总为损伤率估计。

---

## 验证与回归 Validation and Regression

- **回归套件：626 passed / 8 skipped / 0 failed。**
- **1D 冻结锁**：`src/climbing/pad.py` 与其测试保持逐位（bit-identical）不变，`git diff -- src/climbing/pad.py tests/test_climbing.py` 为空。
- **1D/2D 一致性**：默认参数下 2D 轨迹与 1D 逐位相等，2D 扩展不污染已有基线。
- **外部依赖**：载荷链测试依赖 OpenSim 与 FEBio，环境未装时这部分用例无法完整执行。

对标文献的验证靶子（V 编号沿用方案 v2）：

| 靶子 | 内容 | 来源 | 状态 |
|---|---|---|---|
| V1 | 足部骨折临界高度 | [Y25] 群集 2 | 部分达成 |
| V2 | 低于 7 m 无骨折 | [Y25] §2.3 | 达成 |
| V3 | 头部应力峰值时刻 | [Y25] §2.3 | 未达成（绝对值） |
| V4 | 两端先于骨干断裂 | [Y25] Table 1 | 达成 |
| V5 | 损伤部位分布形状 | [B25] Table 1 | 未达成 |
| V6 | 垫刚度与角度-高度转变 | [B25] / [T22] | 达成 |
| V7 | 阶跃 vs 渐进部位分类 | [Y25] §2.3 | 未达成 |
| V8 | 峰值时刻随高度提前 | [T22] Table 3 | 单调性达成，绝对值未达成 |
| V9 | 内翻损伤轻于外翻 | [T22] §结果 | 达成 |

未达成的两类根因：**工况错配**（V3 / V8 绝对值来自刚性地面文献，本模型是软垫）与**模型范围局限**（V5 需要上肢 / 膝 / 肘自由度，本模型没有）。二者都不是调参能解决的，见下节。

---

## 已知局限 Important

这一节是全仓库可信度的落点。模型能做**相对比较与机制分析**，不能做**绝对预测**。

### 1. 绝对骨折高度不可移植

本模型是**单骨截断**模型。同一根骨的「首次骨折高度」会随建模选择（粗体素、平滑网格、真实 AM50 几何、弹性跖面）飘 **3-4 倍**。任何单一数值都不足以作为骨折高度的绝对陈述。

### 2. 上界式 1D 口径的部分量不可引用

1D 通道按上界式口径建立。躯干只有屈曲弹簧一条通路，「腿屈完后躯干直接落在垫上」无法表达，躯干加速度被高估约 **1.5-2 倍**。姿势间的相对排序可用，绝对分档不可引用。

### 3. 工况错配：验证靶子与本模型不同工况

多数验证靶来自**刚性地面**文献（冲击持续 1-5 ms，高度 5-20 m），本项目建的是**软垫**工况（冲击窗 56-67 ms，包线 2-5 m）。因此：

- 绝对峰值时刻结构性不可达（本模型 26.6-37 ms，文献靶 ≤ 12.5 ms）。
- 踝骨折拆分偏差大：模型 22.8% 对实测 8.0%（2.85× 偏高），该拆分被 V6 校准与 `share` 耦合锁死，属结构性问题。

结论：**唯一出路是补软垫物理、做工况对齐，不是继续调参数。**

### 4. 模型范围局限：覆盖不了临床分布形状

模型只覆盖 **7 个骨部位 + 踝 + 头颈**，没有**上肢 / 膝 / 肘**自由度。因此即便软垫物理完全对齐，也无法复现临床的部位分布形状（下肢 67% / 膝 15% / 肘 16%）。蒙特卡洛中 7 个骨部位全为零骨折，踝是唯一非空项，所以「踝第一大」在本模型内属平凡成立，不构成分布证据。

### 5. 有效包线（超出即无效）

| 参数 | 有效范围 |
|---|---|
| 坠落高度 | 2-5 m |
| 落地角 β0 | {10°, 20°, 30°} |
| 滑动摩擦系数 μ_slide | 0.4 |
| ankle_load_share | 0.477 |

**超出包线（如 > 10 m 或 β0 > 30°）外推无效。**

### 可引用与不可引用

| 可引用（带条件） | 不可引用 |
| --- | --- |
| 滑动与卡住相变的存在性，以及角度与高度的转变关系 | 任何处的绝对骨折高度 |
| 内翻轻于外翻的相对关系 | 踝骨折占比 22.8%（实测 8.0%） |
| 峰值时刻随高度单调提前 | 峰值时刻绝对值 26.6-37 ms |
| 踝单项量级 34.5%，乘数 1.16× | 任何 N / MPa 应力绝对值 |
| `rotation='without'` 子集踝 0/321 的结构性对照 | V5 分布形状 |

### 一句话边界

**相对（排序、机制、单项量级）可用；绝对（N、MPa、%、峰值时刻）不可用。**

详细论证与逐项证据见 [`docs/reports/最终结论_2026-10-07.md`](docs/reports/最终结论_2026-10-07.md)。

---

## 文档导览 Documentation

完整索引见 [`docs/README.md`](docs/README.md)。建议按下面顺序阅读：

1. [`docs/reports/最终结论_2026-10-07.md`](docs/reports/最终结论_2026-10-07.md)：当前研究的最终结论与可信边界（先读这份）。
2. [`docs/paper/论文v2_数学建模版.md`](docs/paper/论文v2_数学建模版.md)：论文 v2 稿（数学建模版），模型公式的权威出处。
3. [`docs/paper/论文v2_出处索引.md`](docs/paper/论文v2_出处索引.md)：正文每条数字对应的来源与工况标签。
4. [`docs/methods/方案更新_v2.md`](docs/methods/方案更新_v2.md)：结构性诊断与 P0-P5 重构路线。
5. [`docs/reports/总体进度.md`](docs/reports/总体进度.md)：总体进度总览。

补充材料：`docs/methods/` 收纳各扩展方案与建模交接；`docs/reports/S5_QUOTABLE.md` 是 S5 取数的唯一权威入口；`docs/s5-contact/` 是脚-垫接触收敛系列，成套阅读。

---

## 第三方模型与文献（重要）

`model/`（约 320 MB）与 `paper/`（约 163 MB）**不随本仓库分发**，原因是第三方许可与文献版权。需要自行获取的项：

| 资源 | 来源 | 说明 |
|---|---|---|
| THUMS AM50 V71 Occupant | 丰田（Toyota） | 需申请，附使用条款 |
| OpenSim FullBodyModel | OpenSim 官方 | 见 Rajagopal 等 2016 |
| IITD Ankle-foot model | IITD | 足踝 NLSR 类 `.OSIM` 模型与几何 |
| 文献 PDF | 各出版方 | 按 `docs/paper/` 各篇给出的 DOI / URN 自行获取 |

`temp/`（约 6 GB）是仿真中间产物（`.feb` / `.xplt` / `.hdf5` / `.msh` 等），全部可由脚本再生，无需获取。三者均在 `.gitignore` 中排除，检出后不会自动生成。

---

## 许可与引用 License and Citation

本项目以 **MIT** 许可发布，条款见 [`LICENSE`](LICENSE)。

引用元数据见 [`CITATION.cff`](CITATION.cff)。发布前需注意：其中作者字段为占位（`climbing_fall_analysis contributors`），仓库地址字段也带有 `OWNER` 占位，**填写后才能正式发布**。

可用 `citation-file-format` 工具由 `CITATION.cff` 生成引用：

```bash
cffconvert -f bibtex -i CITATION.cff
```

等价的 BibTeX 形式如下（作者与 URL 待填）：

```bibtex
@software{climbing_fall_analysis,
  title     = {climbing\_fall\_analysis: 室内抱石与攀岩坠落损伤的多层级建模},
  author    = {{climbing\_fall\_analysis contributors}},
  version   = {0.1.0},
  year      = {2026},
  license   = {MIT},
  note      = {Repository URL 与作者见 CITATION.cff（发布前需填写）}
}
```

---

## 贡献与联系 Contributing

欢迎提 Issue 与 Pull Request。报告复现问题或数值偏差时，请附上运行脚本、`pyproject.toml` 中的版本与 `results/` 中对应的输出片段，便于定位。

改动 `src/climbing/pad.py` 前请注意 1D 逐位冻结锁：默认参数下的基线数值不允许漂移，新增能力应走 2D 或独立模块。
