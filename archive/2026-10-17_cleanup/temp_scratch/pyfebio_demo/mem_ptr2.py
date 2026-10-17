"""MEMORY.md 再削：合并踝段 + 纠正韧带线循环论证。"""
from pathlib import Path

M = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md")
s = M.read_text(encoding="utf-8")
REPS = [
    # 合并「刚度2倍」+「F未收敛」+「几何没问题」6 行 → 4 行
    ("""- ⚠️ **刚度 2 倍**：tet4 **259.3** vs penta6 **116.9** N（比值恒 2.2~2.35，五假设已排除，`CONV_TODO.md`）⇒ 疑与接触分辨/罚因子有关
- ⚠️ **F 未收敛但非伪影**：k=F/δ 随罚趋近物理值 **49/151/451/880** N/mm(罚 1e-3→1)，
  对照 Anderson 独立估计 E·A/t≈**1059**（p=1 差 −17%）；但 **p≥2 全负 Jacobian** ⇒ 无法
  验收敛，真值可能高 1.5~2× ⇒ `436.8N` **偏低**。`|u|max` 非判据（悬臂弯 0.17°）
- ✅ **几何没问题**（点在-tet 互穿仅 6 节点）；接触 31mm²@t0→**344mm²**@加载 ≈ Anderson
  生理 300mm²。`[LRN-054/056]`""",
     """- ⚠️ **刚度 2 倍未归因**：tet4 259.3 vs penta6 116.9 N（恒 2.2~2.35，五假设已排除，`CONV_TODO.md`）
- ⚠️ **F 未收敛但非伪影**（`[LRN-054/056]`）：k=F/δ 趋近物理值 **49/151/451/880** N/mm
  (罚 1e-3→1) vs Anderson 独立估计 **1059**；但 **p≥2 全负 Jacobian ⇒ 无法验收敛**，
  真值可能高 1.5~2× ⇒ `436.8N` **偏低**。`|u|max` 非判据（悬臂弯 0.17°）；
  **几何没问题**（互穿仅 6 节点；接触 31→**344mm²**≈生理 300mm²）"""),
    # 纠正韧带线结案
    ("""- **韧带线已降级**：Tochigi 实测站姿韧带本就松弛（CFL 0.2~2.1%）⇒ 我 ΔL≈0 与实验一致；吴恺 1.92~5.24mm 不符实验""",
     """- **韧带线**：原以「ΔL≈0 与 Tochigi 一致」结案 ⇒ 实为**循环论证**（ΔL≈0 由双端夹持强制）⇒ 待重开"""),
]
for a, b in REPS:
    if a in s:
        s = s.replace(a, b, 1)
    else:
        print(f"⚠ 未匹配: {a[:38]}...")
M.write_text(s, encoding="utf-8")
sz = M.stat().st_size
t = M.read_text(encoding="utf-8")
print(f"MEMORY.md -> {sz:,}/10240  余量 {10240-sz:+,}  {'✓ OK' if sz<=10240 else '✗ 超'}")
print(f"  LRN-057: {'Y' if 'LRN-057' in t else 'n'}   循环论证: {'Y' if '循环论证' in t else 'n'}")
