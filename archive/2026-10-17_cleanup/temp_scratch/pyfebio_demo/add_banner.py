"""压缩前收口：给仍含"已作废拉伸侧结论"的文档加顶部更正横幅（纯 prepend，安全）。"""
import io
import os

R = r"D:\Project\climbing_fall_analysis"
BANNER = """> 🚨 **2026-10-05 重大更正（本文件下方部分结论已作废，务必先读这里）**
>
> **1. 加载方向符号一直是反的**（`[ERR-20261005-009]`）：`LOAD_SIGN=-1`（默认）/`DISP_MM>0`
> 是把关节**拉开（拉伸）**，不是压实。**本文件中出现的
> `187.77 / 187.90 / 202.29 / 209.55xx / 592.6x N` 全是拉伸值，不得再作为
> "载荷上限 / 稳定窗 / 600 N 能否达到"的依据。**
> **正解：`LOAD_SIGN=+1` + `DISP_MM<0`（压缩）⇒ `DISP_MM=-1.0` 给 −608.82 N，全程 NORMAL**（余量到 −1017 N）。
>
> **2. "力控制病态 / 关节初始张开（0.97 mm）/ 载荷上限 ~235 N" 三条已作废**：
> 前者是方向符号造成的假象；"1 mm 间隙"是**旧软骨**读数，当前模型**最小间距 0.0211 mm（本就闭合）**；
> 压缩侧可直接到 1017 N。
>
> **3. 权威口径**：`docs/踝FE可用模型使用说明.md`（唯一运行入口）｜
> `[LRN-20261005-072]`（压缩标定曲线）｜`[ERR-20261005-008/009/010/011]`｜`[LRN-20261005-073]`（姿势与判据）。
> **姿势**：倾斜式 `POSE_TILT`/`_ML`（±20°）、几何式 `POSE_IE`（**±5°**，真内外翻）、`POSE_DF`（±2°）。
> **⚠️ `POSE_IE` 打印的"内翻+"是反的：+IE 实际是外翻。**
>
> ---
>
"""
FILES = [
    "ANKLE_FIX_TODO.md", "ANKLE_GOAL_TODO.md", "ANKLE_600N_TODO.md",
    "TIE_PORT_TODO.md", "ANKLE_LIG_TODO.md", "FE_PIPELINE_TODO.md",
    "docs/总体进度.md",
]
MARK = "2026-10-05 重大更正（本文件下方部分结论已作废"
for rel in FILES:
    p = os.path.join(R, rel.replace("/", os.sep))
    if not os.path.exists(p):
        print(f"✗ 不存在: {rel}")
        continue
    t = io.open(p, encoding="utf-8").read()
    if MARK in t:
        print(f"· 已有横幅: {rel}")
        continue
    io.open(p, "w", encoding="utf-8").write(BANNER + t)
    print(f"✓ 加横幅: {rel}  {len(t.encode('utf-8'))} -> "
          f"{len((BANNER+t).encode('utf-8'))} B")
