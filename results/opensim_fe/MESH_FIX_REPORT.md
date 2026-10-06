# FE 网格质量修复报告（gmsh-reparam 干净网格 vs 退化塌缩 tet4）

> 任务：把 4 根 solid-only 骨的退化塌缩 tet4 换成 `part_meshing` 的 gmsh-reparam 干净网格，
> 重跑这 4 行 FE，验证 σ_vm 是否回到物理合理量级。
> **本报告所有数字均为本机实测；失败/发散以 FEBio 原文记录；结论不粉饰。**
> 生成时间：2026-10-04；仓库 `D:\Project\climbing_fall_analysis`；FEBio 4.13.0（`D:\Program\FEBioStudio\bin\febio4.exe`）。

---

## 0. 一句话结论

**干净网格建成且质量/体积都很好，但 σ_vm 的 p95 没有普遍回到物理量级 —— 因为 p95 的主控项是"边界条件伪影"（端部 30% 全约束夹持 + 点节点力 + 斜置骨被沿 bbox 轴斜向加载 → 弯曲），不是塌缩 tet 退化。**

- **跟骨 calcaneus_r**：✅ 成功。σ_vm p95=43.7 MPa，是皮质名义应力 117.4 MPa 的 **0.37×**，首次骨折高度 **1 m → 32 m**。
- **胫骨 tibia_r / 股骨 femur_r**：❌ p95 仍虚高（分别是名义皮质应力的 **31.9× / 18.8×**）。等载荷对比：tibia @5 kN 旧 p95=172 → 新 ≈440–455 MPa（新的反而 **~2.6× 高**，细网格把 BC 奇异点解得更尖）；femur @ref 旧 1169 → 新 708（1.65× 低）。两骨热点都在**中段**（=斜向加载弯曲），机制与旧网格一致——**塌缩 tet 不是 σ_vm 虚高主因；p95 由 BC 决定，网格越细 p95 可能越高**。中位数（体应力）已回到物理量级（0.1–1.6× 名义）。
- **腓骨 fibula_r**：全载荷发散（`negative jacobians`）——细长斜置杆在斜向偏心压缩下**屈曲**（欧拉屈曲应力 ≪ 施加应力），且腓骨实际不独担 25.3 kN。

---

## 1. 做了什么（逐步，命令原文）

| 步 | 动作 | 命令（仓库根，`$env:PYTHONPATH="src"`） | 结果 |
|---|---|---|---|
| 1 | 抽 CORT 外表面 | `extract_part.py --name <b>_cort --pids <CORT_PID> --report-only --out .../parts_extract/<b>_cort` | 4 骨 CORT 面均 `boundary_edges=0` 闭合；tibia CORT 有 16 条 non-manifold 边（2 个部分塌缩 hex） |
| 1b | CORT 面去重修复 | 对 tibia CORT 面按 sorted-face-key 去重（8 个重复面）→ `*_surface_dedup.stl`（0 open / 0 non-manifold） | tibia CORT 可用 raw 路由 |
| 2 | CORT 网格 | `part_meshing.py --char-len {3/4/5}` | 4 骨全部成网格，**体积与原始 CORT 单元体积差 <0.5%**（真·皮质壳） |
| 3a | CORT 网格跑 FE | `run_remesh_fe.py`（单 CORT 实体域，E=18000/ν=0.3） | tibia CORT 因 **3 mm 壁厚在 cl=4 只有 ~1 层单元 → 细长 tet 翻转**，FEBio 原文 `7 negative jacobians detected`；cl≤1.5 mm 才能 >2 层，但会超 5 万 tet 上限 |
| 1' | **回退：抽整骨 union 面**（任务许可） | `extract_part.py --name <b>_union --pids <SPON+CORT> --report-only` | 4 骨 union 面（tibia 17 / femur 2 条 non-manifold → iso-dec；calc/fib raw） |
| 1'b | tibia union 去重 | 去重 10 个重复面 → 0 open/0 non-manifold | tibia union 走 **raw** 路由（min_jac 0.0049→**0.0714**，15× 改善） |
| 2' | union 网格 | `part_meshing.py --char-len {3,4,3,5}` | 4 骨全部 ≤ 2.8 万 tet（≤5 万上限），体积与原始单元体积差 0.45–2.47% |
| 3' | union 网格跑 FE（**主口径，与旧法同 BC**） | `run_remesh_fe.py --bone <b> --load <ref>` | calc/femur @ref 收敛；tibia 降至 10 kN 收敛并线性折回；fibula 全载荷发散 |
| 4 | 物理自检 | σ_nom = load / A（A_cort 与 A_fe 两种） | 见 §3 |
| 5 | 首次骨折高度 before/after | `σ_vm(h)=σ_vm_ref·load(h)/load_ref`，load(h) 直接读 `fracture_matrix.json.loads_n` | 见 §4 |
| 诊断 | PCA 对齐轴向加载（隔离网格 vs BC） | `run_remesh_fe.py --align-pca` | 见 §6 |

**选了哪个表面**：CORT 双壁面能成网格（体积正确），但薄壁 sliver 令 FE 发散；按任务许可改为**整骨（SPON+CORT）并集外表面**成实心皮质骨（E=18000）建模。tibia union 面另做去重以启用 raw 路由。

---

## 2. 干净网格规模 / 质量

### 2.1 实际采用的整骨 union 网格（单实体域，E=18000）

| 骨 | cl (mm) | 方法 | nodes | tets | V_mesh (mm³) | V_raw 单元和 (mm³) | 体积误差 | min_jac | frac jac<0.3 | ≤5万 |
|---|---|---|---|---|---|---|---|---|---|---|
| calcaneus_r | 3 | gmsh-reparam(raw) | 3,664 | 16,845 | 95,642.5 | 96,073.1 | 0.45% | 0.0688 | 2.99% | ✅ |
| tibia_r | 4 | gmsh-reparam(raw, dedup) | 6,513 | 27,508 | 351,562.1 | 354,506.0 | 0.83% | 0.0714 | 3.25% | ✅ |
| fibula_r | 3 | gmsh-reparam(raw) | 3,550 | 11,622 | 50,105.5 | 51,373.0 | 2.47% | 0.0174 | 5.31% | ✅ |
| femur_r | 5 | gmsh-reparam(iso-dec) | 6,101 | 24,653 | 578,837.4 | 586,641.2 | 1.33% | 0.0400 | 3.83% | ✅ |

### 2.2 CORT-only 皮质壳网格（建成但 FE 发散，仅记录）

| 骨 | cl | 方法 | tets | V_mesh | V_cort_raw | 体积误差 | min_jac |
|---|---|---|---|---|---|---|---|
| calcaneus_cort | 3 | gmsh-reparam(raw) | 10,938 | 18,225.8 | 18,279.1 | 0.29% | 0.0140 |
| tibia_cort | 4 | gmsh-reparam(raw, dedup) | 21,710 | 118,551.0 | 118,811.5 | 0.22% | 0.0162 |
| fibula_cort | 3 | gmsh-reparam(raw) | 16,208 | 30,884.2 | 31,015.3 | 0.42% | 0.0141 |
| femur_cort | 5 | gmsh-reparam(raw) | 20,925 | 182,763.3 | 182,990.2 | 0.12% | 0.0249 |

> 网格质量指标：CORT 壳 `min_jac` 0.014–0.025、`frac<0.3` 3.6–7.0%；union 实心 0.017–0.071、`frac<0.3` 3.0–5.3%。均无负体积。
> 注：`<b>_cort_summary.json` 里的 `volume_error_pct`（89%/57%/82%）是 **pyvista `.volume` 对双壁非流形面不可靠**所致；FE 用的 `volume_mm3`（tet 求和）与原始 CORT 单元体积差 <0.5%，已在上表给出。

### 2.3 参考几何（原始 THUMS 单元体积和，用于名义应力）

| 骨 | V_cort (mm³) | L_cort (mm) | **A_cort=V/L (mm²)** | V_union (mm³) | A_union=V/L (mm²) |
|---|---|---|---|---|---|
| calcaneus_r | 18,279.1 | 84.6 | **216.0** | 96,073.1 | 1,135.3 |
| tibia_r | 118,811.5 | 327.6 | **362.6** | 354,506.0 | 1,082.1 |
| fibula_r | 31,015.3 | 304.3 | **101.9** | 51,373.0 | 168.8 |
| femur_r | 182,990.2 | 452.0 | **404.9** | 586,641.2 | 1,297.9 |

> A_cort 来源：**原始 THUMS CORT 单元体积和 ÷ 骨长**（不是发散定理表面体积，后者对非流形面不可靠）。A_fe = 实际 FE 网格体积/长轴跨度。

---

## 3. FE σ_vm 与 σ_nominal 对比（主口径：与旧法相同 BC）

| 骨 | 参考载荷 (N) | 求解载荷 (N) | 折回因子 | σ_vm_ref **p95** | σ_vm_ref median | σ_vm_ref max | σ_nom(A_cort) | σ_nom(A_fe) | **p95/nom_cort** | p95/nom_fe |
|---|---|---|---|---|---|---|---|---|---|---|
| calcaneus_r | 25,354 | 25,354 | 1.00 | **43.65** | 13.09 | 105.65 | 117.38 | 22.31 | **0.37× ✅** | 1.96× |
| tibia_r | 25,266 | 10,000 | 2.53 | **2,224.15** | 38.17 | 9,053.14 | 69.68 | 23.51 | **31.9× ❌** | 94.6× |
| fibula_r | 25,266 | — | — | **FAIL 发散** | — | — | 247.9 | 153.2 | — | — |
| femur_r | 15,228 | 15,228 | 1.00 | **707.60** | 3.86 | 1,492.86 | 37.61 | 11.84 | **18.8× ❌** | 59.8× |

> σ_vm_ref = 收敛载荷下的统计量 × 线性折回因子（线弹性）。tibia 在 1/2/2.5/10 kN 折回得 1984/2294/2363 MPa 有 ±20% 非线性散布（见线性自检）。
> **median（体应力）已回到物理量级**：calc 0.59× / tibia 0.55× / femur 0.10× （相对 A_cort）；fibula 无值。

---

## 4. 首次骨折高度 before / after（p95 口径）

判据 `σ_vm(h) ≥ σ_c`，`σ_vm(h)=σ_vm_ref·load(h)/load_ref`，`load(h)` 逐高度读自 `results/opensim_fe/fracture_matrix.json`（未重跑 OpenSim）。

| 骨 | σ_c (MPa) | before: 旧 p95_ref | before 首次骨折 | after: 新 p95_ref | after 首次骨折 |
|---|---|---|---|---|---|
| 跟骨 calcaneus_r | 150 | 257.15 | **1 m** | 43.65 | **32 m** ✅ |
| 胫骨 tibia_r | 200 | 871.28 | **1 m** | 2,224.15 | **1 m** ❌ |
| 腓骨 fibula_r | 160 | 656.87 | **1 m** | FAIL | **无（发散）** |
| 股骨 femur_r | 220 | 1,168.84 | **1 m** | 707.60 | **1 m** ❌ |

> 只有跟骨首次骨折高度显著后移（1→32 m）。胫/股因 p95 仍被 BC 奇异点主导，仍判 1 m；腓骨无法给出（FE 不收敛）。

---

## 5. 剩余虚高倍数与根因诊断

**剩余 p95 虚高倍数（相对 A_cort 名义）：** calcaneus **0.37×**（已消除）｜ tibia **31.9×**｜ femur **18.8×** ｜ fibula n/a。

**根因（实测证据，非网格塌缩）：**

1. **斜置骨被沿 bbox 最长轴斜向加载 → 弯曲。** PCA 长轴与 bbox 轴的夹角 + 顶/底带质心偏心：

   | 骨 | bbox 跨度 (mm) | bbox 轴 | PCA 长轴 | 与 bbox 轴夹角 | 顶/底带质心偏心 \|ecc\| |
   |---|---|---|---|---|---|
   | calcaneus_r | [70.8, 50.5, 84.2] | Z | [0.402,0.129,0.906] | 25.0° | 56.0 mm |
   | tibia_r | [296.2, 115.1, 327.1] | Z | [0.662,0.132,0.738] | **42.4°** | **311.0 mm** |
   | fibula_r | [275.0, 80.4, 303.6] | Z | [0.658,0.166,0.734] | **42.8°** | **302.6 mm** |
   | femur_r | [450.0, 117.5, 220.0] | X | [0.917,0.084,0.389] | 23.5° | 364.2 mm |

   胫/腓在 X–Z 平面内倾斜 ~42°，把"纵向关节力"沿 bbox Z 施加 = 相对骨轴 42° 斜向 → 巨大弯曲（顶/底带质心相距整个骨长）。

2. **同载荷下新旧 p95 对比（反驳"塌缩 tet 是虚高主因"）：**
   - tibia @5 kN：旧 p95=172.4 MPa（CORT 壳域，粗 THUMS 网格）→ 新 ≈440–455 MPa（干净实心网格）。新的**更高 ~2.6×**，因为更细/更规则的网格把 BC 奇异点解得更尖。
   - femur @ref（15.2 kN）：旧 p95=1168.8 → 新 707.6（1.65× 低）。
   - 两骨热点都在**中段**：旧 tibia @5 kN 的 121 个热点中 117 个在中段（中段均值 436 MPa）；新 tibia @1 kN 的 1376 个热点中 1256 个在中段。→ 弯曲/奇异机制一致，与塌缩 tet 无关。
   - 结论：**塌缩 tet 不是 σ_vm 虚高的主因**；p95 由 BC 决定，网格越细 p95 可能越高（p95 是 BC 奇异敏感指标，不宜作为"网格质量修复"的验证量）。

3. **端部 30% 带全约束 + 点节点力（NodalForce）→ 应力奇异**，`p95` 天然被奇异点主导；干净网格（更规则、更细）把奇异解得更尖，故 fibula/tibia 的 p95 甚至可能**高于**旧的粗 THUMS 网格。
4. **腓骨屈曲**：细长杆（L≈304 mm，A≈169 mm²，r≈3.8 mm）在斜向偏心压缩下欧拉屈曲应力 ≪ 施加应力 → 全载荷发散（`2446 negative jacobians @1 kN`）。物理上腓骨也不独担 25.3 kN。

---

## 6. 轴向加载诊断（PCA 对齐，隔离网格 vs BC）

把 union 网格绕质心旋转，使 PCA 长轴 → +Z 后再 build（load 变为沿骨轴、端带为横截面）：

| 骨 | 求解载荷 (N) | σ_vm p95 | median | max | σ_nom(A_fe) | p95/nom_fe |
|---|---|---|---|---|---|---|
| calcaneus_r | 25,354 | 34.63 | 12.27 | 56.77 | 23.79 | **1.5× ✅** |
| tibia_r | 25,266 | 758.04 | 21.27 | 3,318.21 | 28.98 | 26.2× |
| femur_r | 15,228 | 213.77 | 3.99 | 454.65 | 12.67 | 16.9× |

> 轴向后胫骨 p95 比斜向的 94.6× 明显下降（→26.2×），股骨 59.8×→16.9×；但即使轴向，p95 仍被**端部全约束夹持奇异**抬高。calcaneus 轴向 1.5× 为最好结果。**median 始终物理**（≈0.3–0.7× 名义）。

---

## 7. 失败 / 发散原文（FEBio）

- **tibia CORT 网格**（先试）：`Invalid boundary mesh (overlapping facets) on surface 14 surface 15`（`part_meshing.gmsh_tet_reparam`，iso-dec 后 generate(3)）。
- **tibia CORT FE**：`FEBio ... * ERROR * 7 negative jacobians detected.`（时间步 0.5）。
- **tibia union @25,266 N**：`* ERROR * 7 negative jacobians detected.` → 降载至 10,000 N 收敛。
- **fibula @25,266 N**：`* ERROR * 3 negative jacobians detected.`；@5,000/2,500 N 类似；@1,000 N 网格日志 `* 2446 negative jacobians detected.` / `2446`（细杆屈曲）。
- **fibula 主口径**：全部候选载荷（25,266 / 10,000 / 5,000 / 2,500 / 1,000 N）发散，`fibula_r_remesh_fe.json.status=failed`。

---

## 8. 诚实边界

1. **本修复没有让 σ_vm p95 普遍回到物理量级**。任务假设"塌缩 tet 导致 ~30× 虚高"经实测**不成立为唯一主因**：主因是 **BC（斜向加载 + 端部夹持 + 点载荷奇异）**，旧/新网格同载荷下中段应力基本一致即为证据。只有跟骨（紧凑、倾角小）达标。
2. **p95 是 BC 奇异敏感指标**：网格越干净越细，奇异点解得越尖，p95 可能不降反升。对比"纯净网格质量"应用 **median/体积/体积误差**，而非 p95。
3. **模型是实心皮质骨（E=18000）**，整骨 union 外表面成体，非"皮质壳"；A_fe 是整截面，A_cort 是原始 CORT 体积/长估的皮质承载面积；两种名义都给出，任务主用 A_cort。
4. **腓骨发散有物理成分**（细杆屈曲 + 不独担载荷）；但也含数值脆弱性（同一 1 kN 两次运行结果不一致）。
5. **线性折回误差**：tibia 在 1–10 kN 折回值散布 ±20%（1,984 / 2,294 / 2,363 MPa），非严格线性；旧报告已记录同类 ±8–33%。
6. **未改任何既有模块默认/签名**；未碰 `scripts/ankle_fe/`、`temp/pyfebio_demo/`、`temp/ankle_*`、`FE_PIPELINE_TODO.md`；未覆盖既有 `temp/opensim_fe/fracture_matrix/`、`bone_feb/`、`parts/`。
7. **载荷方向**：旧报告边界 #5 已指出"关节反力是全局纵向，FE 却沿 bbox 最长轴施加"；本报告用 PCA 夹角量化了它（胫/腓 42°），并证明它是 p95 虚高的主因。

---

## 9. 产物

- 本报告：`results/opensim_fe/MESH_FIX_REPORT.md`
- 新产物根目录：`temp/opensim_fe/fracture_matrix_remesh/`
  - `parts_extract/<b>_cort/`：CORT 外表面 STL + summary（含 tibia 的 `*_surface_dedup.stl`）
  - `parts_extract/<b>_union/`：整骨 union 外表面 STL + summary（含 tibia 去重 STL）
  - `meshes/<b>_cort/`、`meshes/<b>/`：干净网格 npz + summary（`gmsh-reparam` 产物）
  - `<b>/<b>_remesh_fe.json`：主口径 FE 结果（含 attempts、线性自检、名义应力）
  - `<b>/<b>_25266N.feb/.log/.xplt` 等：每骨逐载荷 FEBio 输入/日志/结果
  - `<b>_axial/<b>_axial_remesh_fe.json`：PCA 对齐轴向诊断
  - `run_remesh_fe.py`、`hotspot_probe.py`：新驱动/诊断脚本（不改既有模块）

## 10. 复现 & 回归

```powershell
$env:PYTHONPATH="src"
# 1) 抽 CORT / union 外表面（report-only）
& .venv\Scripts\python.exe scripts\opensim_fe\extract_part.py --name tibia_cort --pids 81000700 --report-only --out temp\opensim_fe\fracture_matrix_remesh\parts_extract\tibia_cort
# 2) 干净网格
& .venv\Scripts\python.exe scripts\opensim_fe\part_meshing.py --name tibia_r --stl temp\opensim_fe\fracture_matrix_remesh\parts_extract\tibia_r_union\tibia_r_union_surface_dedup.stl --char-len 4 --out-dir temp\opensim_fe\fracture_matrix_remesh\meshes\tibia_r --ref-volume-mm3 354506.0
# 3) FE（主口径 / 轴向诊断）
& .venv\Scripts\python.exe temp\opensim_fe\fracture_matrix_remesh\run_remesh_fe.py --bone tibia_r --cort-pid 81000700 --load 25266 --a-cort 362.6 --temp-base temp\opensim_fe\fracture_matrix_remesh
& .venv\Scripts\python.exe temp\opensim_fe\fracture_matrix_remesh\run_remesh_fe.py --bone tibia_r --cort-pid 81000700 --load 25266 --a-cort 362.6 --align-pca --temp-base temp\opensim_fe\fracture_matrix_remesh
```

回归（保持 19 passed / 1 skipped 不变）：

```powershell
& .venv\Scripts\python.exe -m pytest tests\test_opensim_fe.py tests\test_febio_run.py -q
# 实测：19 passed, 1 skipped, 1 warning in 95.51s（警告为 .pytest_cache 目录访问权限，与本任务无关）
```
