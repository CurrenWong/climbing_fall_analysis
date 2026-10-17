"""把 Anderson 双论文的结论接进 CONV_TODO / FE_PIPELINE 入口。"""
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
ND = "`docs/Anderson双论文_关键信息.md`"

# ── ① CONV_TODO.md：建议顺序替换为 Anderson 指向 ──
T = ROOT / "CONV_TODO.md"
s = T.read_text(encoding="utf-8")
old = """## 🔜 建议顺序

1. **交换接触面**（最便宜、最可能一击解决）
2. 若还不是 ⇒ 逐单元刚度贡献对比（同位置单元在原模型里单独加载）
3. 若仍查不出 ⇒ **承认第三层**，去找真正的地面真值
   （吴恺模型文件 / 其他足踝 FE 基准数据），
   而不是继续在内部自洽性上打转"""
new = """## 🔜 建议顺序（**已按 Anderson 双论文重排，2026-10-04**）

> ⭐ **真值已到手**：`docs/Anderson双论文_关键信息.md`（Anderson 2006 BMMB + 2007 J Biomech
> 全文重解析）。下面 1–3 全部来自该文已验证做法。

1. **软骨改「沿骨面法向挤出的六面体层」** —— Anderson 全模型 **all-hex**；
   这从根上绕开 tet4 体积锁死 vs penta6 沙漏的二选一（§3.1）
2. **距骨改「弱弹簧约束」+ 解除内翻/内旋约束** —— Anderson 用 500 N/mm（2006）/
   50 N/mm A-P（2007）弱约束距骨，**正式加载时松开内翻/内旋**，
   由胫距关节面自己控制运动学（§2.2）。这才是「外翻+内旋」机理的正确来源
3. **先量关节间隙分布** —— 我可能把软骨做薄了 2 倍（0.703 vs Anderson 1.5~1.7 mm）；
   疑因未摆到负重中和位（§2.1 / §3.3）
4. 若仍差 ⇒ 交换接触面 / 逐单元贡献对比（原方案，降为兜底）

**降级项**：韧带排序复现吴恺 —— Anderson/Tochigi 2005a 用实验证明
「正常踝运动动力学由关节面接触地形主导，**不是韧带**」，
且 Anderson 的验证模型**完全不建韧带**仍达 1.5~15% 精度（§2.3）。"""
assert old in s, "CONV anchor not found"
s = s.replace(old, new, 1)

# 顶部「关键问题」节追加一行真值提示
s = s.replace(
    "**模型的「刚度」有约 2 倍不确定度、且查不出原因 ⇒ 无法判断 600 N 那组结果可不可信。**",
    "**模型的「刚度」有约 2 倍不确定度、且查不出原因 ⇒ 无法判断 600 N 那组结果可不可信。**\n\n"
    "> ⚠️ **2026-10-04 更新**：真值文献已到手（" + ND + "）。\n"
    "> Anderson 的做法（**all-hex 软骨层 + 距骨弱弹簧 + 骨刚性化 + 不建韧带**）\n"
    "> 给出了明确方向 ⇒ **建议顺序已重排**，见本文件末尾。", 1)
T.write_text(s, encoding="utf-8")
print(f"CONV_TODO.md -> {T.stat().st_size:,}")

# ── ② FE_PIPELINE_TODO.md：入口区加一条 ──
FP = ROOT / "FE_PIPELINE_TODO.md"
f = FP.read_text(encoding="utf-8")
if "Anderson双论文" not in f:
    L = f.split("\n")
    hi = next(i for i, l in enumerate(L) if "\U0001f9ed" in l and "接手" in l)
    L.insert(hi + 1, "")
    L.insert(hi + 2,
             "> 🎯 **真值文献（2026-10-04 到手）**：Anderson 2006 BMMB + 2007 J Biomech 全文\n"
             "> → `docs/Anderson双论文_关键信息.md`。要点：**all-hex 软骨层** /\n"
             "> **距骨弱弹簧约束(500、50 N/mm)+ 解除内翻内旋约束** / **骨视为刚体** /\n"
             "> **验证模型不建韧带** / 软骨 E=12 ν=0.42 / 厚度 **1.5~1.7 mm** /\n"
             "> 精度自声明 **10~15%** / 可下载靶子 **Tochigi 2005a 韧带实验**")
    FP.write_text("\n".join(L), encoding="utf-8")
    print(f"FE_PIPELINE_TODO.md -> {FP.stat().st_size:,}（已插指针）")
else:
    print("FE_PIPELINE_TODO.md 已有指针")

# ── ③ 抽查关键结论是否落盘 ──
nd = (ROOT / "docs/Anderson双论文_关键信息.md").read_text(encoding="utf-8")
ck = [("all-hex", "六面体"), ("距骨弱弹簧", "500 N/mm"), ("不建韧带", "neither ligamentous"),
      ("软骨厚度", "1.7 mm"), ("E=12", "12 MPa"), ("Tochigi", "Tochigi"),
      ("50k 收敛", "50,000"), ("精度声明", "10-15%")]
print("\n落盘抽查：")
for name, k in ck:
    print(f"  [{'OK' if k in nd else 'XX'}] {name}")
