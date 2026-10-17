"""修 7 个文件顶部的横幅（改用不可被约定翻转的证据）+ 查其余残留。"""
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"
FILES = ["ANKLE_FIX_TODO.md", "ANKLE_GOAL_TODO.md", "ANKLE_600N_TODO.md",
         "TIE_PORT_TODO.md", "ANKLE_LIG_TODO.md", "FE_PIPELINE_TODO.md",
         "docs/总体进度.md"]

NEW = """> 🚨 **加载方向：本注记写于 2026-10-06（此前竖过两轮相反的旗，见文末"翻转史"）**
>
> **1. `DRIVER=dispx` 下 `DISP_MM > 0` 才是压实（推向关节）**，`DISP_MM < 0` 是**拉开**（拉伸）。
> **依据是三条不可被"约定"翻转的证据**（不是靠反力符号猜的）：
> ① `.feb` 实体：写入 `<value z> = −DISP_MM`；
> ② 材料点几何：距骨在下方（z̄ −709.4 vs 顶面 −653.5），`DISP=+1.0` ⇒ 顶面 dz=**−0.927**（朝距骨）；
> ③ 距骨全固定·无韧带：`+1.0` ⇒ 接触力 **149.59 N**（咬合）；`−1.0` ⇒ **0.0000 N**（分离）。
> **`LOAD_SIGN` 对 `dispx` 分支无效**（实测 z 只由 DISP 的符号决定）。
>
> **2. 因此历史值 `187.77 / 187.90 / 202.29 / 209.5589 / 592.61 N` 是**压实**值（有效），
> 而 `DISP_MM=−1.0 → −608.82 N` 与 `−0.25…−1.6 → −132…−1017 N` 是**拉伸**值（由韧带撑、已作废）。**
> 零回归基准 `DISP_MM=+1.6 ⇒ 592.6020 N` 逐位复现 ✓
>
> **3. 真·压缩曲线**：+0.25/+0.5/+1.0/+1.6/+2.5 mm → 124.0/230.7/418.1/**592.6**/615.4 N
> ⇒ **600 N 在 `DISP_MM≈+1.9 mm` 达成**（"600 N 到不了"作废）。
>
> **4. `POSE_IE` 语义**：`+IE = 内翻`（代码注释正确）。先前"+IE 实际是外翻"已作废
> （`[ERR-20261006-002]`，根因=用数据反推标签的循环论证）。
>
> **5. ⚠️ 翻转史（写下来防止再翻）**
> `[ERR-20261005-009]`：判为"DISP<0 是压实" ⇒ **错**。
> `[ERR-20261006-003]`：改回"DISP>0 是压实" ⇒ **回到最初约定**。
> **⇒ 凡"方向/符号"类结论，必须先给出上述三类证据之一才准写进文档；只凭反力正负或数值大小皆不算。**
>
> ---
>
"""

for rel in FILES:
    p = os.path.join(R, rel.replace("/", os.sep))
    if not os.path.exists(p):
        print(f"✗ {rel}")
        continue
    t = io.open(p, encoding="utf-8").read()
    if "翻转史" in t:
        print(f"· {rel} 已更新")
        continue
    # 找到旧横幅（以 '> 🚨 **2026-10-05 重大更正' 起，到第一个非引用行的 '---' 后）
    m = re.search(r"^> 🚨 \*\*2026-10-05 重大更正.*?^---\n", t, re.S | re.M)
    if m:
        t = t[:m.start()] + NEW + t[m.end():]
        io.open(p, "w", encoding="utf-8").write(t)
        print(f"✓ {rel} 横幅已替换  -> {len(t.encode('utf-8'))} B")
    else:
        print(f"⚠ {rel} 未匹配到旧横幅，改为前置")
        t = NEW + t
        io.open(p, "w", encoding="utf-8").write(t)
        print(f"✓ {rel} 已前置新横幅")

# ── 查其余残留的错误方向论断 ──
print("\n" + "=" * 88)
print("残留扫描（除横幅外，正文里仍声称 DISP<0 是压缩/历史值是拉伸的地方）")
print("=" * 88)
PAT = re.compile(r"(DISP_MM\s*<\s*0|DISP_MM=-(0\.|1)|LOAD_SIGN=\+1|全是拉伸|拉伸值)")
for rel in (FILES + ["docs/踝FE可用模型使用说明.md", "docs/总体进度.md",
                     "ANKLE_FIX_TODO.md", "CONV_TODO.md"]):
    p = os.path.join(R, rel.replace("/", os.sep))
    if not os.path.exists(p):
        continue
    t = io.open(p, encoding="utf-8").read()
    ls = t.splitlines()
    hits = [(i + 1, l.strip()) for i, l in enumerate(ls)
            if PAT.search(l) and i > 25]
    if hits:
        print(f"\n--- {rel}（{len(hits)} 处，正文部分）")
        for i, l in hits[:12]:
            print(f"  {i:>4}: {l[:140]}")
