"""MEMORY 精修至 ≤10240。"""
import io

P = r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md"
t = io.open(P, encoding="utf-8").read()
print(f"改前 {len(t.encode('utf-8'))} B")
REPL = [
    ("- ✅ **A1+A2 已打通**（2026-10-06）：**可用配置** = `SUBTALAR=1` + "
     "**`SUBTALAR_LAUGON=PENALTY`**（AUGLAG 更差）+ `TALUS_SUPPORT=spring` + "
     "`TALUS_K=500` + **`MAX_UPS=0`**（必需；BFGS 开 ⇒ t=0.0045 立崩）。"
     "**可用量程 DISP ≤ +1.2mm（≈519 N）；600 N 够不到**。"
     "@+1.0mm 分解：`MAX_UPS=0` 仅 +0.004%（纯数值工作区）；"
     "**仅接距下 +14.2%｜仅解除约束 −11.7%｜合并 +7.5%**。"
     "⛔ `cart_m1_fix` 反而更差（t 0.76→0.15）**已否证**。"
     "② 力控找根接 OpenSim",
     "- ✅ **A1+A2 已打通**：配置=`SUBTALAR=1`+**`SUBTALAR_LAUGON=PENALTY`**"
     "（AUGLAG 更差）+`TALUS_SUPPORT=spring`+`TALUS_K=500`+**`MAX_UPS=0`**（必需）；"
     "**量程 DISP ≤+1.2mm（≈519 N），600 N 够不到**。@+1.0mm 分解："
     "`MAX_UPS=0` 仅 +0.004%（纯数值）；**仅接距下 +14.2%｜仅解除约束 −11.7%｜合并 +7.5%**"
     "（强交互，不可线性叠加）。⛔ `cart_m1_fix` 更差（t 0.76→0.15）**已否证**。"
     "② 力控找根接 OpenSim"),
    ("- 🔍 **2~4× 的来源**：① **距下关节 = +14.2%**（2026-10-06 实测 @+1.0mm/@1.0 定位移；"
     "旧记的 **44.8% 未获支持、口径不同亦不可直接比**）② 无摩擦切向滑移 **39.5%** "
     "③ tie 罚滑移 **13.9%**；骨 <6%",
     "- 🔍 **2~4× 的来源**：① 距下关节 = **+14.2%**（实测 @+1.0mm；"
     "旧记 44.8% **未获支持**）② 无摩擦切向滑移 **39.5%** ③ tie 罚滑移 **13.9%**；骨 <6%"),
]
for a, b in REPL:
    if a in t:
        t = t.replace(a, b)
    else:
        print("  ⚠ 未命中：", a[:36].replace("\n", "\\n"))
io.open(P, "w", encoding="utf-8").write(t)
n = len(t.encode("utf-8"))
print(f"改后 {n} B  " + ("✅ 余量 %d" % (10240-n) if n <= 10240 else "❌ 仍超 %d" % (n-10240)))
