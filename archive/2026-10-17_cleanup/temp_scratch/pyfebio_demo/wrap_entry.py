"""收口②：在 ANKLE_FIX_TODO.md 顶部插入接手入口。"""
import io
import os

R = r"D:\Project\climbing_fall_analysis"
P = os.path.join(R, "ANKLE_FIX_TODO.md")
t = io.open(P, encoding="utf-8").read()
if "🧭 接手入口（2026-10-06 收口" in t:
    print("· 已有，跳过")
    raise SystemExit(0)

HEAD = """# 🧭 接手入口（2026-10-06 收口 · 压缩后先读这一节）

## 一句话状态
踝 FE 子线**可用于「姿势风险排序」（相对量）**；**绝对数值不可用**（2~4× 不确定度，来源已定位，非未知黑箱）。

## 三件必须知道的事
1. **`dispx` 下 `DISP_MM > 0` = 压实**（三条不可被约定翻转的证据见下方横幅）。
   零回归基准：`+1.6 mm ⇒ 592.6020 N`（逐位复现 ✓）。`[ERR-20261005-009]` 判反、已作废。
2. **相对可用 / 绝对不可用（最重要）**：姿势风险排序经**两套软骨网格交叉验证**
   —— `tibia_cort`/`calcaneus_cort` 排序**完全一致**、`fibula_cort` 首位一致（**+3° 最高**）
   ⇒ **「内翻 ⇒ 腓骨/外踝风险最高」可用**。
   **2~4× 的不确定度来自三条建模伪影（已定量）**：① **缺距下关节接触 44.8%**
   （`subtalar_*` 软骨**已生成但未接进 `.feb`**）② **无摩擦 `sliding-elastic` 接触切向滑移 39.5%**
   （关节间只 0.029 mm 真压扁）③ **tie 罚滑移 13.9%**；**骨合计 <6%、弯曲 <0.2%**。
3. **`[LRN-080]`「韧带承担 77%」已撤回**：实体应变能 147.9 vs 外力功 474.1 N·mm
   ⇒ **~70% 的功在 tie/contact 罚弹簧里** ⇒ 按应变能的分摊不可靠，勿再引用。

## 下一步（按序，含理由）
1. **接入距下关节（subtalar）** —— **最大单一伪影（44.8%）**，而且软骨网格**已经生成**
   （npz 里有 `subtalar_calcaneus` / `subtalar_talus`）⇒ **纯接线、不需新几何**。
   它同时改变**距骨/跟骨的应力路径** ⇒ 是**功能性缺口**，不只是精度问题。
2. **力控找根**（外层割线，压缩曲线单调 ⇒ 3~4 次命中）⇒ 结构上接通 OpenSim
   （现状「OpenSim 给力、模型只吃位移」是断的）。
3. 修软骨厚度上界（max 2.236 mm）／tie 罚因子／接触摩擦。
⚠️ **不要在三条伪影修完前继续雕琢绝对数值**（那是小项上优化）。

## 一键复现（已实测）
```
acceptance.py DRIVER=dispx DISP_MM=+1.6 DET_MIN=0.0 TRUNC_MM=60 AXIS_SRC=full \\
  ALIGN_AXIS=1 FIX_TALUS=transverse CART_NEW=1 TIE=1 CART_DIR=temp/thums/cart_m1 POSE_IE=0
```
⇒ `FEBio rc=0`、**两顶面合力 592.6020 N**、`|u|max = 1.99272 mm`（零回归）。

## 章节索引（本节之后）
| 轮次 | 内容 |
|---|---|
| §第十三 | **柔度归属（三条伪影，subagent 交付）** + 撤回「韧带 77%」 |
| §第十四 | **P1 应力输入独立验证**（重构；tet4 比值 0.98~1.03） |
| §第十五 | subagent `cart_patch.py` v3 解缠验收（≤0: 11→0；但有贴合度回归） |
| §第十六 | **相对结论稳健性验证**（姿势排序跨网格一致） |
| §第七~十二 | P1 骨折裕度首版／P2 软骨网格（三次推翻→定位）／P2'a `area_proj` 修复／P2② 采样诊断 |
| §第一~六 | 符号真相／双靶／C2/C3 收口 |

## 关键文件
`docs/踝FE可用模型使用说明.md`（运行入口）｜`references/ankle-fe-criteria.md`（C2/C3 判据）｜
`docs/Anderson双论文_关键信息.md`｜`CONV_TODO.md`｜
脚本：`scripts/ankle_fe/{ankle_lig_feb,cart_patch,acceptance}.py`｜
诊断：`temp/pyfebio_demo/{p2_locate,p2_cause,p2a_repro,robust_poses,verify_stress3}.py`

---

"""

io.open(P, "w", encoding="utf-8").write(HEAD + t)
n = len(io.open(P, encoding="utf-8").read().encode("utf-8"))
print(f"✓ 已插入接手入口 -> {n:,} B")
