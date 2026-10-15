# AXIAL LOAD AUDIT —— 中轴关节载荷幅值 + `LOAD_SHARE` 来源核对

> **日期**：2026-10-15 ｜ **性质**：只读审计（**零改冻结输入、零改判据、零 FE**）
> **动机**：T6 的 σ 主因在**载荷**不在 BC（`P2_A2_REGFIX_IMPACT.md` §4）——
> σ∝F 严格线性、`A_eff≈93.9 mm²`；要 A1 的窗口打开需把胸椎载荷降 **1.87×**
> **脚本**：`temp/opensim_fe/audit_loadshare.py`、`audit_grf_pad.py`、`temp/find_grf.py`（只读，未提交）

---

## 0. 结论

| # | 结论 | 量级 |
|---|---|---|
| **A** | **`LOAD_SHARE` 与档案记载的来源不符** —— 四个中轴骨的份额**没有一个**对得上 `risk_1d_loadshare` 的 mass-above 比值 | T6 **0.6104 vs 0.5292** |
| **B** | **现役载荷链用「硬地面」而非「软垫」** —— `ground_reaction` 默认 `on_pad=False`（= `hard_surface`，`pad_label="硬地面（离垫）"`），`fracture_matrix.py` 未传 `on_pad` | 场景与论文（抱石+软垫）不一致 |
| **C** | ⚠️ **但换软垫救不了那 1.87×** —— h=5–10 只降 **1.13–1.20×**；**h≥20 软垫反而更差**（垫触底变硬） | 见 §3 表 |
| **D** | 🔴 **真正的量级问题是整条 GRF 标尺高了 10–30 倍** —— 模型 h=1–10 是 **48.8–109.6 ×BW**（硬地面）/ **28.6–91.2 ×BW**（软垫），文献落体实验是 **3.0–4.2 ×BW** | 见 §4 |

**⇒ 「场景用错」不是那个 1.87×。真正指向的是 [D]：有效缓冲行程。**（§5）

---

## 1. 结论 A：`LOAD_SHARE` 来源核对 —— **不符**

`p2_head_contact.py:67` 的 `LOAD_SHARE` 与 `risk_1d_loadshare.py` 的 mass-above 模型对比：

| 骨 | 节段 | 臂 | `LOAD_SHARE` | mass-above 比值 | 差 |
|---|---|---|---:|---:|---:|
| L3 | L3 | 计 | 1.0000 | 0.8954 | +0.1046 |
| **T6** | T6 | 计 | **0.6104** | **0.5292** | **+0.0812** |
| C5 | C5 | 不计 | 0.1840 | 0.1509 | +0.0331 |
| 颅骨 | head | 不计 | 0.1442 | 0.0946 | +0.0496 |

- `risk_1d_loadshare` 的口径：`F_level = (m_above(level)/m_above(L5/S1))·F_lumbar`，
  质量分数取 Pearsall 1996（胸腰椎）+ Ivancic 2006（头/颈）+ Dempster/Winter（上肢 10%），
  参考 `M_REF = 49.70 % BM = 37.44 kg`。
- **档案表述「取自库内 `risk_1d_loadshare`」不准确** —— 四个值全部偏大且**比例不恒定**
  （T6 比 1.153、C5 1.219、颅骨 1.524、L3 1.117）⇒ **不是同一个模型的简单缩放**。

**⇒ 未找到 `0.6104` 的确切出处。当前应在件中把它记为「来源不明（与 `risk_1d_loadshare`
的 mass-above 口径不符）」，而不是「取自 `risk_1d_loadshare`」。**

### 1.1 附：`F_lumbar/GRF` 自洽性检查（**通过，但有一处值得记账**）

刚体上体预期：`F_lumbar/GRF = m_above(L5)/m_total`。THUMS `lumbar` 子树实测质量
**34.24 kg / 75.337 kg = 0.4545**；mass-above 表给 **0.4970**。

| h | GRF 峰值 (N) | F_lumbar (N) | 比值 |
|---:|---:|---:|---:|
| 1 | 32 233.8 | 14 502.7 | **0.450** ✓ |
| 2 | 42 263.7 | 18 761.0 | 0.444 |
| 5 | 53 064.8 | 22 087.6 | 0.416 |
| 10 | 86 149.4 | 32 234.4 | 0.374 |
| 20 | 139 038.0 | 48 810.7 | 0.351 |
| 30 | 180 677.9 | 62 675.7 | 0.347 |
| 50 | 248 330.5 | 86 203.5 | 0.347 |

**⇒ h=1 时 0.450 与刚体预期 0.4545 吻合 ✓**；但**高落高降到 0.347**
⇒ 上体并非刚性减速（可物理解释），**不算错，但应在件中记账**。

---

## 2. 结论 B：现役载荷链用**硬地面**而非**软垫**

```python
# src/climbing/coupling/opensim_grf.py:254
def ground_reaction(*, height_m=5.0, mass_kg=75.337,
                    posture=DEFAULT_POSTURE,     # = "feet-first-stiff"
                    on_pad: bool = False, ...):  # ← 默认 = 刚性地面

# scripts/opensim_fe/fracture_matrix.py (~line 334)
grf = ground_reaction(height_m=float(h), mass_kg=MODEL_MASS_KG)   # ← 未传 on_pad
```

- `pad.py:502`：`on_pad=False` ⇒ `pad = hard_surface(k_n_per_m=1.0e7)`
- `pad.py:723`：`pad_label = "Crash pad" if on_pad else **"硬地面（离垫）"**`
- `pad.py:491` 文档：`True = 落在软垫上；False = **离垫坠落**，改用 hard_surface()`
- `opensim_grf.py` 文档亦标：「`on_pad=False and POSTURES 键=="feet-first-stiff"`」为 L2/L3 门控范围

**⇒ 9 骨骨折矩阵（A1/A2/A3 的输入）跑的是「离垫坠落到硬地面 + 僵硬落地姿势」，
而论文场景是「抱石落地落到软垫上」。** 这是一处**场景/配置不一致**。

⚠️ **但它不是那 1.87×（见 §3）。**

---

## 3. 结论 C（**否掉简单解释**）：换软垫只降 1.13–1.20×，且高落高反而更差

`ground_reaction(on_pad=True/False)` 峰值 GRF 实测（质量 75.337 kg，BW = 739.1 N）：

| h (m) | 硬地面 (N) | ×BW | 软垫 (N) | ×BW | 软垫/硬地面 |
|---:|---:|---:|---:|---:|---:|
| 1 | 36 095.8 | **48.8** | 21 100.6 | **28.6** | 0.58（降 1.71×） |
| 2 | 45 494.0 | 61.6 | 34 248.8 | 46.3 | 0.75（降 1.33×） |
| 5 | 57 074.3 | 77.2 | 50 473.6 | 68.3 | 0.88（降 1.13×） |
| 10 | 80 983.2 | 109.6 | 67 407.2 | 91.2 | 0.83（降 1.20×） |
| 20 | 141 816.0 | 191.9 | 173 282.7 | 234.5 | **1.22（更高）** |
| 30 | 184 026.3 | 249.0 | 261 858.0 | 354.3 | **1.42（更高）** |
| 50 | 251 632.1 | 340.5 | 395 032.9 | 534.5 | **1.57（更高）** |

**⇒ 在真正关心的 h=5–10，换软垫只给 1.13–1.20×，远不够 1.87×；h≥20 软垫反而更差**
（垫厚 0.2 m 在高能量下触底 ⇒ 底刚段 `bottoming_stiffness=50` 起作用）。

**⇒ 「现役用错场景」成立（结论 B），但不能解释 T6 的 σ 与 A1 的窗口缺口。**

---

## 4. 结论 D（**真问题**）：整条 GRF 标尺高了 10–30 倍

| | 落高 | 质量 | 峰值 vGRF | ×BW |
|---|---:|---:|---:|---:|
| **文献**（同行评审落体实验，屈腿柔性落地） | 0.3 m | 65 kg | 1 910 N | **3.0** |
| **文献** | 0.5 m | 65 kg | 2 655 N | **4.2** |
| **文献** | 0.5 m | 70 kg | 2 710 N | **3.9** |
| 现役模型（硬地面） | 1 m | 75.3 kg | 36 096 N | **48.8** |
| 现役模型（软垫） | 1 m | 75.3 kg | 21 101 N | **28.6** |

**⇒ 落差 10–30 倍，且随落高持续放大。** 无论走硬地面还是软垫路径，**都远超文献量级**。

---

## 5. ★ 机制定位：有效缓冲行程（§4 的成因）

能量守恒（三角脉冲近似）：`mgh ≈ ½·F_peak·s` ⇒ `s ≈ 2mgh / F_peak`

| 场景 | mgh | F_peak | 推算行程 `s` |
|---|---:|---:|---:|
| 文献 0.5 m / 65 kg | 319 J | 2 655 N | **240 mm**（≈24 cm 屈膝） |
| 现役 h=1 m | 736 J | 36 096 N | **41 mm** |
| 现役（仅硬面项 `k_n=1e7 N/m`） | 736 J | — | **3.6 mm** |

**⇒ 文献的柔性落地有 ~24 cm 缓冲行程；现役模型只有 4–41 mm ⇒ 有效刚度高约 10 倍。**

**指向的具体位置（未验证）**：

- `POSTURES` 默认 = **`feet-first-stiff`**（僵硬落地）。可选项里另有
  `controlled-drop` / `toe-point` / `tuck-roll` 等更柔的姿势。
- `simulate_boulder_fall` 的 `k_flex`（人体屈曲刚度）—— 是行程的主控量。
- `hard_surface` 的 `k_n_per_m=1.0e7` —— `pad.py` 自认「**全模型最弱参数**」。

**⇒ 下一手应当是：① 核对 `feet-first-stiff` 是不是抱石软垫场景的正确姿势；
② 审计 `k_flex` 与 `k_n` 的出处；③ 用「缓冲行程」做独立校验（文献 ~24 cm）。**

---

## 6. 诚实边界

1. **零 FE、零改冻结输入、零改判据**；全部产出为只读审计与新件。
2. **§4 的文献对比有口径差异**：文献是**受控柔性屈腿落地**（被试主动缓冲），
   而现役模型的 `feet-first-stiff` 是**僵硬落地**。**两者不完全可比** ——
   但 10–30 倍的落差**不是姿势差异能解释的量级**（僵硬落地文献上界约 10 ×BW）。
   **建议**：补一次**僵硬落地**条件下的文献/实测锚点再下定量结论（**未做**）。
3. **§3 的 `on_pad=True/False` 对比是同一模型内的相对比较**，未与文献核对绝对值。
4. **`LOAD_SHARE = 0.6104` 的确切出处未找到**（§1）；本件只证「与档案所述来源不符」。
5. **未改任何既有判定**：A1 仍 FAIL、A2 仍 4/5（现役）。本件**不产生**新的判定。
6. **§5 的行程反推是量级估算**（三角脉冲假设），非实测。

---

## 7. 溯源

| 项 | 位置 |
|---|---|
| `LOAD_SHARE` | `scripts/opensim_fe/p2_head_contact.py:67` |
| mass-above 模型 | `scripts/opensim_fe/risk_1d_loadshare.py`（`SEG_PCT` / `MASS_ABOVE` / `RATIO`） |
| 静止载荷链调用 | `scripts/opensim_fe/fracture_matrix.py`（`ground_reaction` 调用处 ~L334） |
| `ground_reaction` 签名与默认 | `src/climbing/coupling/opensim_grf.py:254-266` |
| `on_pad=False` ⇒ hard_surface | `src/climbing/pad.py:502`；标签 `:723`；说明 `:491` |
| 姿势表 | `src/climbing/pad.py:260`（`POSTURES`） |
| 载荷链缓存（GRF/lumbar） | `temp/opensim_fe/fe_bc_probe/load_chain.npz` |
| 文献锚点（落体峰值 vGRF） | CliffsNotes 生物力学实验数据表（0.3/0.5 m, 65/70 kg） |
| 上游动机 | `results/opensim_fe/P2_A2_REGFIX_IMPACT.md` §4 |
| 审计脚本 | `temp/opensim_fe/audit_loadshare.py`、`audit_grf_pad.py`（只读，未提交） |

**复现**：

```powershell
cd D:\Project\climbing_fall_analysis
& .venv\Scripts\python.exe temp\opensim_fe\audit_loadshare.py    # LOAD_SHARE 来源 + F_lumbar/GRF
& .venv\Scripts\python.exe temp\opensim_fe\audit_grf_pad.py      # 硬地面 vs 软垫峰值 GRF
```
