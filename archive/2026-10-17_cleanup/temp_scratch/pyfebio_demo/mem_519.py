"""MEMORY：把 A1+A2 那条压缩，腾出空间写入「519 N 是收敛代价非承载极限」的机制。"""
import io

P = r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md"
t = io.open(P, encoding="utf-8").read()
n0 = len(t.encode("utf-8"))
print(f"改前 {n0} B")

OLD = ("- ✅ **A1+A2 已打通**：配置=`SUBTALAR=1`+**`SUBTALAR_LAUGON=PENALTY`**"
       "（AUGLAG 更差）+`TALUS_SUPPORT=spring`+`TALUS_K=500`+**`MAX_UPS=0`**（必需）；"
       "**量程 DISP ≤+1.2mm（≈519 N），600 N 够不到**。@+1.0mm 分解："
       "`MAX_UPS=0` 仅 +0.004%（纯数值）；**仅接距下 +14.2%｜仅解除约束 −11.7%｜合并 +7.5%**"
       "（强交互，不可线性叠加）。⛔ `cart_m1_fix` 更差（t 0.76→0.15）**已否证**。"
       "② 力控找根接 OpenSim")
NEW = ("- ✅ **A1+A2 已打通**：配置=`SUBTALAR=1`+**`SUBTALAR_LAUGON=PENALTY`**"
       "（AUGLAG 更差）+`TALUS_SUPPORT=spring`+`TALUS_K=500`+**`MAX_UPS=0`**；"
       "**量程 ≤+1.2mm（519 N），600 N 够不到**。@+1.0mm 分解：`MAX_UPS=0` 仅 +0.004%"
       "（纯数值）；**仅接距下 +14.2%｜仅解除约束 −11.7%｜合并 +7.5%**（强交互）。"
       "⛔ `cart_m1_fix`/`cart_m3` 均已否证。② 力控找根接 OpenSim\n"
       "- 🧭 **「519 N」是收敛代价、不是承载极限**（2026-10-07 同位移实测）：位移驱动 ⇒ "
       "519 N 只是最后一个能解出的点；**A2 解除约束后距骨转动 0.52°（基线 0.0033°，158×）"
       "⇒ 胫距软骨由均匀受压变边缘受压 ⇒ 薄楔单元(1.2mm)被剪翻**。"
       "对照：基线(锁死)可达 **615.4N@+2.5mm** ⇒ 挪向物理正确，收敛余量掉 ~16%。"
       "⚠️ 勿写成「踝只能承 500N」")
if OLD in t:
    t = t.replace(OLD, NEW)
    print("✓ 已替换 A1+A2 条并加入机制")
else:
    print("⚠ 未命中")
io.open(P, "w", encoding="utf-8").write(t)
n1 = len(t.encode("utf-8"))
print(f"改后 {n1} B  " + ("✅ 余量 %d" % (10240-n1) if n1 <= 10240 else f"❌ 超 {n1-10240}"))
