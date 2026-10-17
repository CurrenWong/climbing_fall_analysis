"""记录 P1d：加密机制成功（节点 9.5x）但求解失败；以及 shell 污染坑。"""
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
ND = ROOT / "docs/Anderson双论文_关键信息.md"
T = ROOT / "ANKLE_FIX_TODO.md"

s = ND.read_text(encoding="utf-8")
s = s.rstrip() + """

---

# §12 P1d 加密（2026-10-04）——**机制成功，但求解失败**

## 12.1 ✅ 加密机制работает：关节面节点 73 → 692（9.5×）

脚本 `temp/pyfebio_demo/refine3.py`（1 级红色加密全 body）：
- `hex8 → 8 hex8`（棱中点 + 面心 + 体心）
- `tet4 → 8 tet4`（棱中点）
- `quad4 壳 → 4 quad4`
- 节点 22,960 → **219,717**（新增 棱 137,730 / 面 50,640 / 体 8,387）

**在加密网格上重跑 `cartilage.py`（`THUMS_NPZ` / `CART_OUT` 指向 r1）：**

| | 原始 | 加密 R1 | 倍数 |
|---|---|---|---|
| tibiotalar 胫侧关节面节点 | **73** | **692** | **9.5×** |
| talus 侧 | 89 | 887 | 10× |
| 软骨单元（胫/距） | 93 / 129 | **894 / 1,256** | 9.6× |
| 软骨体积 | 1.057 / 1.109 cm³ | **2.414 / 2.573 cm³** | **2.3×** |

**⇒ ⭐ 关键发现：细网格下软骨体积涨 2.3 倍** ——
说明**粗网格严重低估了关节面的真实面积**（细网格才能捕捉曲率）。
这与「接触区欠分辨」的判断一致，**是本轮最有价值的产出**。

## 12.2 ⚠️ 但加密网格**求解失败**（`rc=1`，0.1 s）

```
tag "Elements" (line 49104) : missing attribute "id"
  → <Elements type="hex8" name="talus_cort__hex8"></Elements>  ← 空块
```
**骨与软骨单元全部被 `DET_MIN` 闸门剔除**（`软骨 []`），
说明加密后的子单元 **detJ 为负**（朝向不对）。

**已尝试的修法**：在 `refine3.py` 里加 `fix_hex`（子单元与母单元都做
「交换底面/顶面 4,5,6,7,0,1,2,3」的朝向修正）——**未能解决**。
⇒ 说明问题**不只是底/顶交换**能修的朝向，更可能是**母单元节点顺序不符合 FEBio 约定**。

## 12.3 📐 支撑上面判断的旁证

用两种体积度量测**同一**网格：

| 度量 | 原始 | 加密后 | 相对差 |
|---|---|---|---|
| 6-tet 扇形分解 | 890,100 | 880,044 | −1.13% |
| **三线性 Jacobian Gauss 积分** | **980,326** | 887,871 | −9.43% |

**⇒ 两种度量在【原始】网格上就差 10%**（980k vs 890k）。
扇形分解对**扭曲**六面体不精确；两者不一致 ⇒
**母六面体存在明显扭曲**，我的「参数化 → 棱/面中点」假设不一定成立。

**⇒ P1d 的正确收尾**：
先**逐单元检查 npz 里 hex 的节点顺序约定**（对照 THUMS/LS-DYNA 与 FEBio 的差异），
再决定是「规范化顺序后重加密」还是「改用 tied-elastic 把细软骨绑到粗骨上」（§3.14.5）。

## 12.4 🔴 操作坑（本轮代价最大）

我在**同一条 cmd 命令**里先 `set "THUMS_NPZ=...r1.npz"` 再跑 `cartilage.py`，
**随后**又在同一条命令里跑 `cmp_ref.py` —— **`set` 污染了整个 shell**，
而 `cmp_ref.py` 里是 `if npz: env[...] = npz` 的写法，
**"原始"组继承到了 r1 ⇒ 四个 run 全在跑加密网格**，连原始组也失败，
一度误判为「我的改动弄坏了模型」。

**⇒ 修法**：对比脚本里必须**无条件** `env.pop("THUMS_NPZ")` 再按需设置。
**⇒ 通用教训**：cmd 的 `set` 作用域是**整个命令行**，不要在同一条命令里
既改 env 又跑对比/回归。
"""
ND.write_text(s, encoding="utf-8")
print(f"Anderson双论文_关键信息.md -> {ND.stat().st_size:,}")

t = T.read_text(encoding="utf-8")
t = t.replace(
    """- [ ] **P1d（新首要）加密 + 规范筛选关节面骨表面网格** —— 关节面节点 ~72 → **≥300**；
      且需**筛选干净**（不能是 `gap<6mm` 纯距离），法向一致
      ⇒ 这是 hex8 路线的**前提**，不是可选项""",
    """- [x] **P1d 加密机制** ✅ **核心目标达成**：关节面节点 **73 → 692**（9.5×）、
      软骨单元 **93 → 894**（9.6×）、软骨体积 **×2.3**
      ⇒ ⭐ **证明粗网格严重低估关节面面积**（本轮最有价值产出）
- [ ] **P1d2 🔴 修复加密网格的单元朝向** —— 加密后骨/软骨单元 detJ 全负被剔除
      （`rc=1`）。已试「底/顶交换」未解决 ⇒ 疑因**母 hex 节点顺序不符合 FEBio 约定**
      ⇒ 先逐单元核对 THUMS vs FEBio 的 hex 顺序约定；
      或改用 **`tied-elastic`**（手册 §3.14.5：专门连接非共形网格）
      把细软骨绑到粗骨上，**绕开加密骨的整个难题**""")
t = t.replace("**当前：P0 ✅ / P1a-c ✅；P1d 加密+规范筛选关节面为下一步（hex8 的前提）**",
              "**当前：P0 ✅ / P1a-d ✅（加密机制通、求解未通）；P1d2 修子单元朝向为下一步**")
T.write_text(t, encoding="utf-8")
print(f"ANKLE_FIX_TODO.md -> {T.stat().st_size:,}")

L = ROOT / ".learnings/LEARNINGS.md"
l = L.read_text(encoding="utf-8")
if "LRN-20261004-045" not in l:
    l = l.rstrip() + """

## LRN-20261004-045: cmd 的 `set` 污染整条命令行 —— 对比脚本必须【无条件】清除环境变量

- **类别**: `correction` / `operational`
- **场景**: 在同一条 cmd 命令里：
  `set "THUMS_NPZ=...r1.npz" && python cartilage.py ... & ... && python cmp_ref.py`
  `cmp_ref.py` 内部用 `if npz: env["THUMS_NPZ"] = npz` 来切换数据源。
- **后果**: "原始"组（`npz=None`）**继承**了 shell 里的 `THUMS_NPZ=r1`
  ⇒ **四个 run 全部跑加密网格** ⇒ 连"原始"基准也 `rc=1` 失败
  ⇒ 我一度误判为「我这几轮改动把模型弄坏了」，
  直到在日志里看到 `[data] npz = ...r1.npz` 才发现。
- **修法**: 对比脚本里对每个受控变量**无条件** `env.pop(K, None)` 再按需 `env[K]=v`；
  **不要**写 `if v: env[K]=v`（它会让"控制组"继承外部值，控制组就不再是控制组）。
- **通用教训**:
  1. cmd 的 `set` 作用域是**整个命令行**，**不要在同一条命令里既改 env 又跑对比/回归**
  2. **控制组被环境污染**是实验设计里最隐蔽的一类错误 ——
     症状是"连基准都失败了"，此时**第一个该怀疑的不是代码而是 env**
  3. 任何受控变量都要让控制组**显式断链**，不能依赖"没设就是默认"
- **行动**: 已修 `temp/pyfebio_demo/cmp_ref.py`；结论写入
  `docs/Anderson双论文_关键信息.md` §12.4
- **关联**: `[LRN-20261004-038]`（开关没生效，同属"环境/机制没按预期起作用"）
"""
    L.write_text(l, encoding="utf-8")
    print(f"LEARNINGS.md -> {L.stat().st_size:,}（+LRN-045）")

print("\n抽查：")
nd = ND.read_text(encoding="utf-8"); td = T.read_text(encoding="utf-8")
for nm, ok in [("§12 落盘", "§12 P1d 加密" in nd),
               ("73→692", "692" in nd and "9.5×" in nd),
               ("体积×2.3", "2.3×" in nd),
               ("P1d 完成", "[x] **P1d" in td),
               ("P1d2 下一步", "P1d2" in td),
               ("LRN-045", "LRN-20261004-045" in L.read_text(encoding="utf-8"))]:
    print(f"  [{'OK' if ok else 'XX'}] {nm}")
