"""落盘 A1+A2 结果：可用配置、分解、边界、根因、以及对我错误结论的最终替换。"""
import io
import os

R = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"

# ---------- 1) MEMORY：替换那两条（44.8% 与「下一步」） ----------
MP = os.path.join(W, "MEMORY.md")
m = io.open(MP, encoding="utf-8").read()
m0 = len(m.encode("utf-8"))

OLD1 = ("- 🔍 **2~4× 的来源（⚠️ 44.8% 那条未定）**：① **缺距下关节 44.8% —— 未证实也未否证**"
        "（2026-10-06 的「+0.0000%」测量**因 case 未显式传 `SUBTALAR`、实际没接上而作废**；"
        "接入后能否收敛**尚未解决**）② 无摩擦切向滑移 **39.5%** ③ tie 罚滑移 **13.9%**；骨 <6%")
NEW1 = ("- 🔍 **2~4× 的来源**：① **距下关节 = +14.2%**（2026-10-06 实测 @+1.0mm/@1.0 定位移；"
        "旧记的 **44.8% 未获支持、口径不同亦不可直接比**）② 无摩擦切向滑移 **39.5%** "
        "③ tie 罚滑移 **13.9%**；骨 <6%")
if OLD1 in m:
    m = m.replace(OLD1, NEW1); print("✓ MEMORY 伪影① 已替换")
else:
    print("⚠ MEMORY 伪影① 未命中")

OLD2 = ("- **⏭️ 下一步**：① **A1+A2 已跑 4 格矩阵**（E0/E1/E2/E3）⇒ 基线 592.6020 ✓，"
        "E1/E2/E3 **均无负 J、但在 t=0.475/0.585/0.95 收敛停滞（残差卡 ~1e2）**；"
        "`MAX_UPS=0` 已消掉 t=0 失控 ⇒ **下一步是迭代/步长**（`TIME_STEPS`↑、`MAX_REFS`↑）"
        "② 力控找根接 OpenSim")
NEW2 = ("- ✅ **A1+A2 已打通**（2026-10-06）：**可用配置** = `SUBTALAR=1` + "
        "**`SUBTALAR_LAUGON=PENALTY`**（AUGLAG 更差）+ `TALUS_SUPPORT=spring` + "
        "`TALUS_K=500` + **`MAX_UPS=0`**（必需；BFGS 开 ⇒ t=0.0045 立崩）。"
        "**可用量程 DISP ≤ +1.2mm（≈519 N）；600 N 够不到**。"
        "@+1.0mm 分解：`MAX_UPS=0` 仅 +0.004%（纯数值工作区）；"
        "**仅接距下 +14.2%｜仅解除约束 −11.7%｜合并 +7.5%**。"
        "⛔ `cart_m1_fix` 反而更差（t 0.76→0.15）**已否证**。"
        "② 力控找根接 OpenSim")
if OLD2 in m:
    m = m.replace(OLD2, NEW2); print("✓ MEMORY 下一步 已替换")
else:
    print("⚠ MEMORY 下一步 未命中")
io.open(MP, "w", encoding="utf-8").write(m)
n1 = len(m.encode("utf-8"))
print(f"MEMORY {m0} -> {n1} B " + ("✅" if n1 <= 10240 else f"❌超 {n1-10240}"))

# ---------- 2) ROADMAP ----------
RP = os.path.join(R, "ANKLE_ROADMAP_TODO.md")
ADD = """

---

# ✅ A1+A2 已打通（2026-10-06 深夜，主 agent 亲测）

## 1. 可用配置（这是本轮最重要的产物）

```
SUBTALAR=1              接距下关节
SUBTALAR_LAUGON=PENALTY  ★ 关键：AUGLAG 更差（0.46→0.65 的差别）
TALUS_SUPPORT=spring     ★ A2：横向硬约束换成弱弹簧
TALUS_K=500              Anderson 值
MAX_UPS=0                ★ 必需：BFGS 开 ⇒ t=0.0045 立刻崩
CART_DIR=temp/thums/cart_m1（不要用 cart_m1_fix，见 §4）
DISP_MM ≤ +1.2mm
```

## 2. 收敛边界（DISP=dispx，POSE_IE=0）

| DISP_MM | F | 结果 |
|---|---|---|
| +0.25 | 139.98 | ✅ |
| +0.5 | 248.95 | ✅ |
| +1.0 | 449.63 | ✅ |
| +1.1 | 484.87 | ✅ |
| **+1.2** | **518.88** | **✅ 上界** |
| +1.3 / +1.6 | — | ❌ 负 Jacobian |

⇒ **可用量程 ≈ +1.2 mm ⇒ ~519 N；而 600 N 需要更大位移 ⇒ 够不到。**
（基线可达 +1.6mm/592.6N；接 A1+A2 后量程收窄。）

## 3. ⭐ 正确分解（DISP=+1.0 mm，五格全 NORMAL）

| case | 配置 | F | 相对基线 |
|---|---|---|---|
| D0 | 基线（BFGS 开） | 418.097 | — |
| D0n | 基线 + `MAX_UPS=0` | 418.113 | **+0.004%** |
| **D1** | **仅接距下关节（A1）** | **477.387** | **+14.18%** |
| **D2** | **仅解除横向约束（A2）** | **369.029** | **−11.74%** |
| **D12** | **A1+A2 合并** | **449.631** | **+7.54%** |

**结论**：
1. **`MAX_UPS=0` 无物理影响（+0.004%）** ⇒ 可安全当数值工作区用
2. **距下关节贡献 +14.2% 刚度**（正效应，符合「多一条承托路径」）
3. **解除横向约束 −11.7%**（负效应，硬约束→软弹簧）
4. **两者有强交互**：合并 +7.54% ≫ 朴素相加 +2.44% ⇒ **不能线性叠加**

## 4. ⛔ `cart_m1_fix` 已否证（不要用它救收敛）

| 配置 | 结果 |
|---|---|
| 基线 + fix 网格 | rc=0 ✅ 但 **F=656.05 vs 592.60 ⇒ +10.7%** |
| A1 + fix 网格 | ❌ t=0.273（对照 cart_m1 的 0.645） |
| A1+A2 + fix 网格 | ❌ t=0.153（对照 cart_m1 的 0.761） |

⇒ **修翻转单元的网格反而让收敛更早失败**（与它已知的「贴合度回归」一致：
接触更早介入）。**它也不适合当基线（力偏 +10.7%）。**

## 5. 根因链（本轮查清）

| 层 | 现象 | 结论 |
|---|---|---|
| ① 求解器 | `SUBTALAR=1` + BFGS 开 | **t=0.0045 立崩（9 负 J）** ⇒ 拟牛顿更新病态 ⇒ **必须 `MAX_UPS=0`** |
| ② 接触 | subtalar 用 AUGLAG | 卡在 `augmentation # 3`，只到 t=0.463 ⇒ **`PENALTY` 更好（t=0.645）** |
| ③ 单元 | +1.6mm 时 | **单个单元 8241**（`tibiotalar_talus__penta6`，gauss pt 5）翻负 ⇒ 量程上界 |

**③ 的意义**：卡点**不在距下关节本身**，而是**加了新载荷路径后，
原有的胫距软骨（基线里本来就能收敛的那部分）被压过头**。
⇒ 若要突破 +1.2mm，方向是**胫距软骨的网格/厚度/材料**，不是继续调接触参数。

## 6. 我前几轮的错误结论（已就地替换，不叠加）

- ~~「接入距下关节效应 = +0.0000%」「44.8% 被否证」~~ ⇒ **实验错标**
  （`a1_matrix2.py` 的 case 没传 `SUBTALAR`，默认已是 0 ⇒ 那个 case 根本没接上）
  ⇒ 真实值 **+14.18%**
- ~~「`MAX_UPS=0` 后双侧接入完全收敛」~~ ⇒ 同一错标；真实情况是**还需 PENALTY + spring，
  且量程只到 +1.2mm**

## 7. 新增的检查器教训

- 我的负 J 检查器写的是 `(\\d+) negative jacobians?`（复数带计数），
  而 FEBio 实际消息是 **`Negative jacobian was detected at element N at gauss point M`**（单数无计数）
  ⇒ **漏检**，一度让我误判「无负 J」。**已改用两种写法都匹配。**
- **教训**：`-output_negative_jacobians=N` 才是权威读数，别用日志字符串猜。

## 8. ⏭️ 下一步（更新）

1. **突破 +1.2mm 上界** ⇒ 针对 **`tibiotalar_talus` 软骨**（不是距下关节）：
   厚度/网格/材料，或接触参数针对该对调优
2. **力控找根**（B1）：现在有了可用配置，可在外层做二分 ⇒ 逼近 519 N 以内的目标力
3. **接 OpenSim**（D）：需先解决 600 N 够不到的问题（否则场景受限）
"""
if "A1+A2 已打通（2026-10-06 深夜，主 agent 亲测）" not in io.open(RP, encoding="utf-8").read():
    io.open(RP, "a", encoding="utf-8").write(ADD)
    print(f"ROADMAP +{len(ADD)} 字符 -> "
          f"{len(io.open(RP, encoding='utf-8').read().encode('utf-8')):,} B")
else:
    print("· ROADMAP 已有")

# ---------- 3) 日记 ----------
DP = os.path.join(W, "memory", "2026-10-06.md")
DADD = """
## 🎉 A1+A2 打通（2026-10-06 深夜）

**可用配置**：`SUBTALAR=1` + **`SUBTALAR_LAUGON=PENALTY`** + `TALUS_SUPPORT=spring` +
`TALUS_K=500` + **`MAX_UPS=0`**，`CART_DIR=cart_m1`，**DISP ≤ +1.2mm**。

**边界**：+1.2mm ⇒ 518.88 N 上界；+1.3 起失败。**600 N 够不到**（基线可到 +1.6mm/592.6N）。

**正确分解 @+1.0mm（五格全 NORMAL）**：
- `MAX_UPS=0` 自身 **+0.004%**（纯数值）
- **仅接距下关节 +14.18%** ｜ **仅解除横向约束 −11.74%** ｜ **合并 +7.54%**（强交互）

**根因链**：① BFGS 病态（必须 `MAX_UPS=0`）② subtalar 用 AUGLAG 卡在 augmentation ⇒ 改 PENALTY
③ 上界由**单个单元 8241**（`tibiotalar_talus__penta6`）翻负决定 ⇒ **瓶颈是胫距软骨，不是距下关节**。

**⛔ `cart_m1_fix` 已否证**（收敛更早失败 t 0.76→0.15；且力偏 +10.7%）。

**⚠️ 我前几轮的错误（最终替换）**：「接入 +0.0000%」与「44.8% 否证」都源于
`a1_matrix2.py` 的 case **没传 `SUBTALAR`**（当时默认已是 0）⇒ 那个 case 根本没接上。
真实值 **+14.18%**。**教训：跑矩阵前先确认每个 case 的 env 真的生效**
（用「输出里有没有 `[subtalar] 接入`」这类可观测证据自证，而不是信变量名）。
"""
if "A1+A2 打通（2026-10-06 深夜）" not in io.open(DP, encoding="utf-8").read():
    io.open(DP, "a", encoding="utf-8").write(DADD)
    print(f"日记 +{len(DADD)} 字符")
else:
    print("· 日记已有")
