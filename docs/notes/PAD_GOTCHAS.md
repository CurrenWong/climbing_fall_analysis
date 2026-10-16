# `pad.py` / 1D-2D 链模型 · GOTCHAS

> 与 `FE_BC_GOTCHAS.md`（FE/BC 专用）同类，但范围是 **`src/climbing/pad.py`、`pad2d.py`、
> `chain2d*.py`、`bodies.py`** 这条落地动力学线。
> 新增条目请附**实测数字**与**复现脚本**。

---

## 1. 🔴 `peak_force_n`（`F_pad`）对 `k_flex` **不敏感** —— 不得当 GRF / 腿刚度观测量

**现象（硬地面 `on_pad=False`，@0.42 m、80 kg，`anchored-*` 三档）**：

| 档 | `k_flex` (kN/m) | **`F_pad` (kN)** | `F_flex` (kN) |
|---|---:|---:|---:|
| soft | 6.3 | **8.61** | 2.27 |
| self | 8.2 | **8.61** | 2.34 |
| stiff | 18.5 | **8.63** | 2.78 |

**⇒ `k_flex` 变 2.9 倍，`F_pad` 只动 0.2% ⇒ 它由 `hard_surface(k_n=1e7)` × `m_low`
的**刚性接触尖峰**决定，不是腿刚度的响应量。**

- ✅ **腿刚度的观测量 = `result.meta["f_flex_peak_n"]`**（对 k 敏感：2.27→2.78，+22%）。
- ⛔ 用 `peak_force_n` 做的「模型 vs 力台 GRF」对比**会得到无意义的结论**
  （2026-10-16 曾据此误判「模型 @0.42 m 高 3.3×、缺屈服吸能」——已作废）。
- 复现：`scripts/fp_anchor_probe.py`、`tests/test_climbing.py::TestAnchoredPostures::
  test_f_pad_insensitive_to_k_flex`。

---

## 2. 🔴 **L2 刚性门按姿势的「键」判定 ⇒ 传 `Posture` 对象会静默绕门**

`pad.py` 的 `RIGID_STANDING_CHAIN = {"kf_mult": 2.0, "x_max": 0.015, "k_stop": None}`
（2026-10-08 L2 标定采纳）**只在 `on_pad=False` 且姿势的*键* == `"feet-first-stiff"`** 时生效
⇒ 把 `k_flex` **×2.0**（1.1e6 → 2.2e6）、行程限到 **15 mm**。

实现（`simulate_boulder_fall`）：

```python
post = POSTURES[posture] if isinstance(posture, str) else posture
# L2 门用 POSTURES 的**键**判姿势（post.name 是中文显示名）
posture_key = posture if isinstance(posture, str) else None
```

**⇒ 传 `replace(POSTURES["feet-first-stiff"], ...)` 这类**对象**时 `posture_key=None`
⇒ 门不生效 ⇒ 同名姿势两种调用结果不同。**

**实测（`feet-first-stiff` @1 m，80 kg，硬地面）**：

| 调用 | `F_pad` (kN) | `F_flex` (kN) |
|---|---:|---:|
| `posture="feet-first-stiff"`（**键**，含门） | **37.46** | **31.10** |
| `posture=<Posture 对象>`（**无门**） | **33.65** | **27.43** |

**⇒ 差 13.4%。** 后果与用法：

1. **引用既有数字必须标「含门/无门」**：`V1_WHY_20KN_AND_GAP_LOCUS.md` §3 表里标着
   「1.1e6（模型 `feet-first-stiff`）」的 27.43 kN 是**无门**值；真实 `feet-first-stiff`
   在刚性口径下是 **2.2e6 / 15 mm**。
2. **新姿势的键若不是 `"feet-first-stiff"` ⇒ 门不生效**（设计如此：「pad 路径与其他姿势
   逐位不变」）。**照抄姿势字段造新档时，务必确认你要的是有门还是无门口径**
   （`anchored-*` 就是要**无门**，因为比的是 `k_flex` 本身 —— 预注册 §4-O5）。
3. 守卫测试：`tests/test_climbing.py::TestAnchoredPostures::test_l2_rigid_gate_key_vs_object_differs`
   （**特征化**测试：锁定这不是 bug、锁死 ratio 区间）。

---

## 3. 🟡 `POSTURES` 的键 ≠ `Posture.name`

`Posture.name` 是**中文显示名**（如「脚先落但腿绷直」），**不是**键。
任何按名字分支的逻辑都应在**键**上做（同 §2）；日志/报告里出现中文名时不要反查成键。

---

## 4. 🟡 `x_flex_max` 不是主控参数，行程由 `x = F/k` 决定

`x_flex_max` 只是**行程上限**。在触到上限之前，实际行程 `x = F/k` 由刚度决定
（例：`k=1.1e6` ⇒ `x≈25 mm` ≪ 30 mm 上限 ⇒ 上限根本没被触及）。
**⇒ 对 `x_flex_max` 单变量做灵敏度往往会得到"无影响"的假阴性**
（2026-10-16 实测：`feet-first-stiff` 基座 0.03→0.20 m，`F_flex@4 m` 45.98→60.30 kN **不降反升**）。

---

## 5. 🔴 `hard_surface` 的 `k_n` **不是**总刚度 —— 有效刚度 = `k_n × A(κ)`

`CrashPad.force_n(compression, area) = stress_pa(xi) × area`，而
`hard_surface(k_n_per_m=...)` 把该形参喂给 **`linear_modulus_pa`（Pa）**、`thickness=1`
⇒ **有效总刚度**

```
k_eff = k_n × A(κ)          # A(κ) = Posture.contact_schedule 的接触面积
```

**实测（`anchored-stiff` @0.42 m、80 kg）**：

| `k_n` | `F_pad` (kN) | `x` (mm) | **实测 `k_eff` = F/x** | `k_eff/k_n` |
|---|---:|---:|---:|---:|
| 1.0e7 | 8.625 | 15.728 | **548.4 kN/m** | 0.055 |
| 3.0e6 | 5.194 | 29.333 | 177.1 kN/m | 0.059 |
| 1.0e6 | 3.870 | 57.236 | 67.6 kN/m | 0.068 |

**⇒ `hard_surface` docstring 声称的「默认 1.0e7 N/m（动态人体等效刚度）」比**实际**硬 ~18×。**

**用法纪律**：
1. **引用硬地面接触硬度时必须报实测 `k_eff`（= `F_pad/max_compression_m`），不得引用 `k_n` 本身。**
2. 想调接触硬度就调 `RIGID_CONTACT["k_n_per_m"]`（**L3 门**，默认 `None` = 1.0e7）——
   ⚠️ **`k_eff` 是 `k_n` 的非线性函数**（`A(κ)` 随压缩变化，`k_eff/k_n` 在 0.055–0.16 间漂移），
   所以**不能**按 `k_eff` 直接线性反解；须实测。
3. 特征化/守卫测试：`tests/test_climbing.py::TestRigidContactGate`。
4. 锚点与实测见 `results/FP_CONTACT_REPORT.md`、`docs/methods/P3_足部接触柔度_预注册.md`。

---

**最近更新**：2026-10-16（`LEG_K_ANCHOR` → `FP_CONTACT`）。
**锚点**：`src/climbing/pad.py`（`RIGID_STANDING_CHAIN` / `RIGID_CONTACT` / `hard_surface` / `_anchored_posture`）、
`scripts/fp_anchor_probe.py`、`scripts/fp_contact_probe.py`、`tests/test_climbing.py`、
`results/FP_ANCHOR_REPORT.md`、`results/FP_CONTACT_REPORT.md`
