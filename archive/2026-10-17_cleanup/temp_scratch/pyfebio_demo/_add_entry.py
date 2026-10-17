"""FE_PIPELINE_TODO.md -> 顶部加自包含「接手入口」，并修正过时的状态行。"""
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis\FE_PIPELINE_TODO.md")
t = P.read_text(encoding="utf-8")
before = len(t.encode("utf-8"))

old_head = """# FE_PIPELINE — THUMS 骨材料 → 软骨生成 → FEBio 冒烟（按 c→b→a）

用户要求按 **c → b → a** 顺序执行。
**状态：C ✅ 完成 / B ✅ 完成 / A ⚠️ 进行中（.feb 语法全通，卡在 Jacobian）**

---
"""

new_head = """# FE_PIPELINE — THUMS 骨材料 → 软骨生成 → FEBio 求解（c→b→a）

**状态：C ✅ / B ✅ / A ✅ 已收敛出物理结果**（2026-10-04）

---

# 🧭 接手入口（先读这一节就够，其余是细节）

## 一句话

**踝关节（胫距）有限元模型已跑通 FEBio 求解并产出可信应力**；对标吴恺 2012 的
**量级验证已通过**。剩余工作都是**加料**（多关节、韧带），不是**修地基**。

## 三件必须知道的事

| | 内容 |
|---|---|
| ✅ **已验证** | `DRIVER=disp DISP_MM=0.1472` → `rc=0` / `NORMAL TERMINATION` / 16.6 s，支反力 **436.8 N**；软骨 von Mises p99 **2.58 / 4.27 MPa** → 外推 600 N 胫侧 ≈ **3.4 MPa**，与吴恺的 **3.74 MPa** 吻合 |
| ❌ **负面结果（别再试）** | 「提高软骨网格质量」**不是**收敛瓶颈。四道质量闸门让四面体 \\|V\\|min 提升 **42 倍**（0.0013→0.053 mm³）、穿透仍 0，**但收敛更差**（0.30 mm/60 步 0.965→0.642；工作点 17 s→60 s）。⇒ 瓶颈在**接触**。`CART_QUALITY=1` 才开，默认关 |
| 🚧 **已知阻塞** | subtalar 关节（`FULL_ANKLE=1`）**网格建得好**（35,095 单元、骨单元 0 剔除、2 个接触界面注册成功）但**第一个载荷步失败** —— 接触面**非合法 2-流形**（`surf_subtalar_calcaneus` 787 面，边使用次数 `{1:333, 2:592, 3:163, 4:76, 5:9, 6:1}`，合法曲面应全 =2；开放边界 ≈ 被丢弃的 66 个单元挖的洞） |

## 一键复现（已验证配置）

```cmd
cd /d D:\\Project\\climbing_fall_analysis
.venv\\Scripts\\python.exe scripts\\ankle_fe\\thums_materials.py   :: C 材料
.venv\\Scripts\\python.exe scripts\\ankle_fe\\cartilage.py         :: B 软骨
set "DRIVER=disp" && set "DISP_MM=0.1472" && set "TIME_STEPS=20" && ^
  .venv\\Scripts\\python.exe scripts\\ankle_fe\\build_feb.py       :: A 求解（约 20 s）
.venv\\Scripts\\python.exe temp\\pyfebio_demo\\a6_reaction.py     :: 读支反力 + 应力
```

> ⚠️ cmd 的 `set X=Y &&` **行尾不能带空格**（否则 40013，见 MEMORY gotcha #5）。
> `DRIVER=force`（600 N 面力）**会失败** —— 胫骨只靠接触约束 ⇒ 刚体模态 ⇒ 奇异。
> 现在用 `DRIVER=disp` + 顶面全约束（夹持端），**不依赖接触来稳定**。

## 下一步（按价值排序，含理由）

1. **加 8 组踝韧带** ⭐ **建议先做这个**
   数据已在 `temp/thums/thums_lowerlimb.npz`（`81100101~81101401_shell`）。
   理由：**不经过接触这条路径** ⇒ 不会被下面的 2 卡住；
   而且能复现吴恺的**定性排序**（**胫跟＞胫舟＞跟腓＞下胫腓前**），
   这比 3.74 MPa 的绝对值**验证价值更高**（项目计划里验证目标就是「量级 + 定性排序」，
   量级已拿到，排序还是空白）
2. **修 subtalar 接触面 2-流形**
   要点：**先构造完整偏移曲面（合法 2-流形），再做单元筛选**；或对接触面做边界修补。
   验收用 `temp/pyfebio_demo/diag_facets.py`（数每条边被几个面用）
3. **换回 force 驱动（真 600 N 面力）**
   需先解决刚体模态：弱弹簧、或初始 tied 接触

## 环境开关（全部有默认值，直接跑即可）

| 开关 | 默认 | 文件 | 作用 |
|---|---|---|---|
| `DRIVER` | `force` ⚠️ | `build_feb.py` | `disp`=顶面全约束+指定位移（**能收敛**）；`force`=面力（会失败）|
| `DISP_MM` | `0.05` | `build_feb.py` | `disp` 模式的压缩量；**已验证 0.1472** |
| `TIME_STEPS` | `20` | `build_feb.py` | 载荷步数；接触问题细化有效 |
| `FULL_ANKLE` | `0` | `build_feb.py` | 1=接入距下关节（多接触对 + 固定面改跟骨底）|
| `LAUGON` | `PENALTY` | `build_feb.py` | `AUGLAG` 试过，**第一个步就压翻 4 个，更糟** |
| `NODE_RELOC` / `TWO_PASS` / `PENALTY` | `0`/`0`/`1` | `build_feb.py` | 接触参数 |
| `CART_QUALITY` | `0` | `cartilage.py` | 1=开网格质量闸门（**结论是负面的**，见上）|
| `CAP_SAFETY` | `0.97` | `cartilage.py` | 每侧软骨厚度被对侧 offset 面限厚，防初始过盈（**有效，默认开**）|
| `BONES_ONLY` | `0` | `build_feb.py` | 二分定位用（只跑骨、无软骨无接触）|

## 关键文件

| 文件 | 作用 |
|---|---|
| `scripts/ankle_fe/thums_extract.py` | 流式解析 219 MB LS-DYNA `.k` → `thums_lowerlimb.npz` |
| `scripts/ankle_fe/thums_verify.py` | 提取验证（单位/体积/共形/装配） |
| `scripts/ankle_fe/thums_materials.py` | C：解材料卡 → `bone_materials.json` |
| `scripts/ankle_fe/cartilage.py` | B：骨面法向 offset 生成 4 层软骨 |
| `scripts/ankle_fe/build_feb.py` | A：组装网格 + 写 `.feb` + 调用 FEBio |
| `temp/pyfebio_demo/a6_reaction.py` | 读 `.xplt` → 支反力 + 各域 von Mises |
| `temp/pyfebio_demo/diag_facets.py` | **判定接触面是否合法 2-流形**（数边使用次数）|
| `temp/pyfebio_demo/diag_contact.py` | 核对接触面法向朝向 |
| `docs/踝关节FE建模交接.md` | 全程历史（§1–§13，46 KB）|
| `.learnings/LEARNINGS.md` | `LRN-20261003-016~020` = 本轮全部教训 |
| `references/fe-modeling-gotchas.md` | FE 陷阱速查（工作区）|

## 别再踩的坑（全部已实测，细则见 §A）

1. **FEBio `.feb` 语法**：`SolidDomain` 的 `mat` 填**材料名**（不是数字 id）；
   `<Surface>` 子标签**按单元类型命名**（`<quad4>`/`<tri3>`）；
   **段序** `Mesh→MeshDomains→LoadData→Loads→Boundary→Contact→Output→Control`
   —— 顺序错会被报成 `unrecognized tag` 这种**误导性错误**。完整模板见本文件末尾
2. **只靠接触约束的部件 ⇒ 刚体模态 ⇒ 刚度矩阵奇异** ⇒ 位移爆掉、单元成批翻负
   （实测 2200 个）。→ **把加载端整个截面全约束**（夹持端）
3. **tet4 的"翻转"必须是奇置换** —— 4 节点整体倒序 = 3 次交换 = **偶置换** ⇒
   det 符号**不变**！正解：**交换两个节点**。曾因此误剔 122/462 + 183/618 个软骨单元
4. **接触面必须与单元网格的真实面重合** —— 否则 FEBio 接触**完全不检测**，
   但求解报 `rc=0`/`NORMAL TERMINATION`、**无任何警告**，**只是应力全为 0**（极误导）
5. **接触面还必须是合法 2-流形** —— 每条边恰好被 2 个面用
6. **初始过盈让接触"从第一步激活"会直接压翻单元**（t=0 就有接触力）
7. **二分定位别用零载荷变体** —— 残差本已 8e-14，但 FEBio 判据要求降到 0，
   永远满足不了 ⇒ 报 `Problem is diverging` ⇒ **假失败**
8. **别输出 `<var type="contact traction"/>`** —— 写出的 `PLT_FACE_DATA` 会让
   `pyfebio.xplt.to_hdf5` 崩掉；且它对自定义 Surface 集合会 `KeyError`（需容错补丁）

---
"""

assert t.count(old_head) == 1, f"锚点命中 {t.count(old_head)}"
t = t.replace(old_head, new_head)
P.write_text(t, encoding="utf-8")
print(f"FE_PIPELINE_TODO.md {before:,} -> {len(t.encode('utf-8')):,} bytes")
