# I″ · `share = 0.2241` 真解 FE —— **FAIL**，I′ 的落窗无法确认

> **日期**：2026-10-13 ｜ **驱动**：`scripts/opensim_fe/i3_fibula_real_solve.py`（纯增量，**未改** G9 预注册常量）
> **产物**：`results/opensim_fe/g9_fibula_fe_support_share02241.json`、`G9_FIBULA_FE_SUPPORT_REPORT_share02241.md`
> **FEBio 日志**：`temp/opensim_fe/g9_fibula_fe_support/share_0.2241__long95/fib_g9.log`

---

## 0. 一句话

**`share = 0.2241`（5662 N）在 `long95` 侧向支撑带下 FEBio 发散（238 negative jacobians）⇒ I′ 的 `OR = 1.2701` 落窗无法确认；
但同件的 0.10 档复现成功 ⇒ 腓骨线本身可复现。**

---

## 1. 两个结果

### 1.1 `share = 0.10` 同条件复现 —— **PASS**（附带的重要正面结果）

| | G9 记录（`..._d1a2.json`） | 本件真解 |
|---|---|---|
| `σ_nom/σ_c` 下端 | 10.53 MPa | **10.5 MPa** |
| `σ_nom/σ_c` 上端 | 76.93 MPa | **76.86 MPa** |

**⇒ 逐值吻合（~0.3%）。腓骨线可复现。**
**⇒ 与冻结 BC2 线形成对照**：那条线同条件重跑差 **6–22%**（`FE_BC_PROBE.md`）。
**⇒ 「不可复现」是**冻结 BC2 线**的特有缺陷，不是腓骨 FE 的通病。**

### 1.2 `share = 0.2241`（5662 N）真解 —— **FAIL**

```
share=0.2241__long95: FAIL (build_run) BoneFebError: FEBio 失败（退出码 1）
* 238 negative jacobians detected.
* Problem is diverging. Stiffness matrix will now be reformed
------- failed to converge at time : 1 / 0.833333 / 0.666667 / 0.5 / 0.333333 / 0.166667
```

**⇒ 从 `time = 0.166667`（即载荷的 1/6）起，每一个时间步都发散，负 Jacobian 稳定在 238–239 个。**
**⇒ 不是收敛器参数问题、不是偶发**：`long95` 已把屈曲压掉（0.10 档线性 ratio 2.106），
但载荷上到 5662 N 时进入**大变形/网格失稳**。

---

## 2. ★ 失败签名与本项目三处已知失败同族

| 处 | 签名 |
|---|---|
| **本件** `share=0.2241` | **238–239 negative jacobians**，逐步发散 |
| 踝 tied 罚线 | 「**罚 ≥ 2 全负 Jacobian**」 |
| F-1（股骨） | 侧向屈曲型不稳定 ⇒ ≥ 某载荷发散 |
| T6（胸椎） | 参考载荷 22.1 kN 即发散 |

**⇒ 四者是同一类现象：刚性/大变形约束下的接触面或侧向失稳，而非各自独立的随机失败。
⇒ 这也是本项目制度「同源归因」的一次正面应用。**

---

## 3. ⇒ 对 I′ 结论的影响（**如实回退，不保留 PASS 候选**）

| | I′ 的说法 | 本件之后 |
|---|---|---|
| `share=0.2241 ⇒ OR=1.2701 落窗` | 「PASS 候选」 | **降级为【未验证的外推值】** —— 外推区间 2.24× 的**上端点无法解出**，该值既未被证实也未被证伪 |
| **A3 判定** | 有望 3/3 | **维持 2/3**（不变） |

**⇒ 腓骨仍未命中，但原因现在更清楚：不是「share 取不到」，
而是「能解出的 share 上限（≤ 0.10 档、2527 N）对应的 OR（1.1229）恰好差窗口下沿 4.6%，
而能落窗的 share（≥ 0.2241）恰好解不出来」。**
**⇒ 二者之间隔着一条数值可解性边界。**

---

## 4. 诚实边界

1. **verdict 块打印的 `None` 是取不到主设定造成的假阴性**：
   G9 的 `PRIMARY = ("share_0.050", "mid")` 指向一个本件未跑的组合
   （我把 `SHARES` 缩到两档）。**`None` 不代表解算失败**，真实状态见 §1。
2. **本件未改 G9 的 `SHARES` 常量**，是运行时追加；产物走 `--out-suffix _share02241`，**未覆盖**既有预注册产物。
3. **未试其它支撑带**（`long85` / `three`）或更高 `share` 的中间档（如 0.15 / 0.18）——
   **「可解上限」目前只知道落在 (0.10, 0.2241] 之间，未定位。**
4. **未改动任何已签字判定**；A3 现状 **2/3**、A2 现状 **阶跃 5/6 + 渐进 0/3 = 5/9** 均不变。

---

## 5. 下一步（二选一，待拍板）

| # | 项 | 说明 |
|---|---|---|
| **a** | **二分定位可解上限** | 在 (0.10, 0.2241] 间二分（如 0.15 / 0.18）⇒ 给出腓骨可解 share 上限与该处 OR。**零新增设计、纯取证** |
| **b** | **攻可解性** | 参照 F-1 的做法（侧向覆盖 ≥80% 自由长 + 夹持 DOF 无关）再加强支撑，或改加载方式（跟垫式分布载荷）⇒ 把 0.2241 推出来 |

**⇒ 建议先做 (a)**：它能**判定「落窗口是否落在可解区内」**——
若可解上限对应的 OR < 窗口下沿，则 (b) 才有必要；
若可解上限就能落窗，则腓骨直接 3/3，不需要 (b)。

---

## 6. 溯源

| 项 | 来源 |
|---|---|
| 本件驱动（纯增量） | `scripts/opensim_fe/i3_fibula_real_solve.py` |
| 被复用的 G9 逻辑（零改动） | `scripts/opensim_fe/g9_fibula_fe_support.py`（`SHARES` 常量原样，仅运行时追加） |
| 本件产物 | `results/opensim_fe/g9_fibula_fe_support_share02241.json`、`G9_FIBULA_FE_SUPPORT_REPORT_share02241.md` |
| FEBio 失败日志 | `temp/opensim_fe/g9_fibula_fe_support/share_0.2241__long95/fib_g9.log` |
| 前两件 | `I_FIBULA_SHARE_AXIS.md`（冻结行）、`I2_FIBULA_SHARE_AXIS_SUPPORTED.md`（侧向支撑行） |
| 同族失败签名 | 踝 tied 罚线（MEMORY）、`FE_BC_STAB.md`、`T6_BC_DIAGNOSIS.md` |