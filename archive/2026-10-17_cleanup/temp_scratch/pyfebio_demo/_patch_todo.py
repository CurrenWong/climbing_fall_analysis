from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis\FE_PIPELINE_TODO.md")
t = P.read_text(encoding="utf-8")

old = """## A 的下一步（可选优化）
- 提高软骨网格质量（厚度更均匀，消除 det(J)<0.1 的薄片）→ 应能改善收敛包络
- 换 force 驱动（600 N 面力）：需先解决刚体模态，可用弱弹簧或 `tied` 初始接触
- 加 subtalar 关节（现已生成 4 层软骨，但当前模型只用了 tibiotalar 两层）"""

new = """## A 的下一步 —— 第二轮结果（2026-10-04）

- [x] 提高软骨网格质量 —— **试过了，假设是错的**。质量闸门（`CART_QUALITY=1`）把
      四面体 `|V| min` 从 0.0013 提到 0.053 mm³（**42 倍**）、穿透仍 0，
      **但 FEBio 收敛反而变差**（0.30 mm/60 步：0.965 → 0.642；工作点 17 s → 60 s）。
      ⇒ **收敛瓶颈是接触，不是薄片单元**。见 `LRN-20261003-020`
- [x] 接入 subtalar 关节 —— **网格建好了，接触面不合法（未完成）**。
      `FULL_ANKLE=1` 可建出 35,095 单元 / 10,324 节点（骨单元 0 剔除）、
      两个接触界面都注册成功，但**第一个载荷步就 `2 negative jacobians`**（载荷≈0，
      细化载荷步无效）。
      根因：接触面**不是合法 2-流形** —— `surf_subtalar_calcaneus` 787 个三角面，
      边使用次数 `{1:333, 2:592, 3:163, 4:76, 5:9, 6:1}`（合法曲面应全部 =2）。
      开放边界数 ≈ 被丢弃的 66 个软骨单元挖出的洞。FEBio 也警告 `18 invalid facets`。
      **注意**：单关节（能收敛）模型同样有 2/3 个 invalid facet 却能跑完
      ⇒ 警告本身不致命，致命的是非流形的**规模**。见 `LRN-20261003-019`
      → **下一步**：先构造完整偏移曲面（合法 2-流形）**再**筛单元，
      或对接触面做边界修补/重网格
- [ ] 换 force 驱动（600 N 面力）：需先解决刚体模态（弱弹簧或初始 tied 接触）
- [ ] 加 8 组踝韧带（`81100101~81101401_shell` 已在 npz 里）→ 能复现吴恺的
      **定性排序**（胫跟＞胫舟＞跟腓＞下胫腓前），比 3.74 MPa 绝对值更有验证价值

### 第二轮新增的机制
| 开关 | 位置 | 作用 |
|---|---|---|
| `CART_QUALITY=1` | `cartilage.py` | 打开软骨网格质量闸门（默认关，结论见上） |
| `CAP_SAFETY=0.97` | `cartilage.py` | 每侧软骨厚度再被**对侧已生成的 offset 面**限制，防初始过盈（默认开） |
| `FULL_ANKLE=1` | `build_feb.py` | 接入距下关节（跟骨 + 第二对接触 + 固定面改到跟骨底） |
| `LAUGON=AUGLAG` / `NODE_RELOC=1` / `TWO_PASS=1` / `PENALTY=x` | `build_feb.py` | 接触算法与参数（试过 AUGLAG+node_reloc，第一个步就压翻 4 个，更糟） |
| `TIME_STEPS=n` | `build_feb.py` | 载荷步数（细化对首步冲击有效，但对"载荷≈0 就失败"的问题无效） |

### 诊断脚本
- `temp/pyfebio_demo/diag_facets.py` —— 解析 `.feb` 的 Surface，统计
  **每条边被几个面使用**（判定 2-流形）+ 面积/重合/重复面。**这是本轮最关键的工具**
- `temp/pyfebio_demo/diag_contact.py` —— 核对接触面法向朝向
- `temp/pyfebio_demo/a6_reaction.py` —— 读支反力（标定载荷量级）

### 单关节模型回归（改动后仍成立，2026-10-04）
`DRIVER=disp DISP_MM=0.1472 TIME_STEPS=20` → **rc=0 / converged / 16.6 s**，
支反力 **436.8 N**（改前 440.7 N），软骨 p99 **2.58 / 4.27 MPa**（改前 2.48 / 4.42），
差异 <1% ⇒ 上面这些改动**没有破坏已验证结果**。"""

assert t.count(old) == 1, f"anchor x{t.count(old)}"
P.write_text(t.replace(old, new), encoding="utf-8")
print("FE_PIPELINE_TODO.md 已更新")
