"""落盘 ERR-20261005-010（_outer 新格式 KeyError + 旧 .feb 陷阱）并查看 MEMORY 尺寸。"""
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"
ERR = os.path.join(R, ".learnings", "ERRORS.md")
SNIP = """
---

## [ERR-20261005-010] `_outer()` 只认旧软骨键 ⇒ `POSE_TILT` 一直崩；且被"旧 `.feb` 陷阱"掩盖

- **症状 A**：任何带 `POSE_TILT` / `POSE_TILT_ML` 的运行 `rc=1`，
  `KeyError: 'tibiotalar_tibia_offsets is not a file in the archive'`。
- **根因 A**：`_outer()`（姿势的关节面间距自检用）只读**旧格式**键 `{j}_offsets`；
  而 `CART_NEW`（`cart_patch.py`）用 `{j}_nodes` + `{j}_is_offset`（节点键名 `cn:{j}:{k}`）。
  ⇒ **姿势功能自新软骨格式起就是坏的**（`[LRN-20261005-063]` 之后的所有姿势结论都不可靠）。
  修法：`CART_NEW` 分支走 `cn:{j}:{k}` + `CN_ISOFF[j][k]`。
- **症状 B（更隐蔽）**：8 个不同 `POSE_TILT` 跑出**完全相同**的 `Fz = −608.82 N`、剪切恒 0。
- **根因 B**：脚本没查**构建返回码**、也没**删旧 `.feb`**。构建失败后
  **FEBio 会读上一次留下的 `.feb`** ⇒ 所有姿势其实跑的是同一个文件。
- **教训（两条都要写进 runner 模板）**：
  1. **构建完必须查 `rc`，并确认 `.feb` 的 mtime 是本次的**；
  2. **每次跑前删旧 `.feb` 与旧读数** —— "结果逐字相同"永远是第一嫌疑。
- **另一条同源教训**：**反力必须看符号**（压缩为负）。只看幅值会把**拉伸**当成成功
  （`[ERR-20261005-009]`）。
"""

txt = io.open(ERR, encoding="utf-8").read()
if re.search(re.escape("ERR-20261005-010"), txt):
    print("★ 已存在，跳过")
else:
    io.open(ERR, "a", encoding="utf-8").write(SNIP)
    print(f"ERRORS.md {len(txt)} -> {len(txt)+len(SNIP)}")

M = r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md"
s = io.open(M, encoding="utf-8").read()
print(f"\nMEMORY.md 字节数 = {len(s.encode('utf-8'))} (上限 10240)")
i = s.find("踝关节 FE 子线")
print("踝段起始位置:", i)
print(s[i - 900:i + 400] if i > 0 else "(未找到)")
