"""补最后 75 B：紧凑化 Anderson 真值行 + 姿势行。"""
import io

MEM = r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md"
m = io.open(MEM, encoding="utf-8").read()
n0 = len(m.encode("utf-8"))
REPS = [
    ("- **Anderson 真值**（`[LRN-067]` 已核实原文）：all-hex｜距骨**弱弹簧 500 N/mm 接地**+解除内外翻内旋｜",
     "- **Anderson 真值**（`[LRN-067]` 已核原文）：all-hex｜距骨**弱弹簧 500 N/mm 接地**+解除内外翻内旋｜"),
    ("  **骨刚性 + reference node（⇒ 才有转动 DOF 可解除；可变形模型做不到）**｜不建韧带｜E=12 ν=0.42、厚 1.5~1.7mm；",
     "  **骨刚性 + reference node**（⇒ 才有转动 DOF 可解除）｜不建韧带｜E=12 ν=0.42、厚 1.5~1.7mm；"),
    ("- ⛔ **「不用软骨」≠解法**（`[LRN-062]`）：力变**罚因子伪影**。\n"
     "  **通用判据：力随罚因子显著变化 ⇒ 不可信**",
     "- ⛔ **「不用软骨」≠解法**（`[LRN-062]`）：力变**罚因子伪影** ⇒ "
     "**通用判据：力随罚因子显著变化即不可信**"),
    ("  （剪切 0→79.7 N）+ **几何式 `POSE_IE` ±5°（真内外翻）/`POSE_DF` ±2°**",
     "  （剪切 0→79.7 N）+ **几何式 `POSE_IE` ±5°/`POSE_DF` ±2°**"),
]
for old, new in REPS:
    if old in m:
        m = m.replace(old, new)
    else:
        print(f"⚠ 未命中：{old[:38]}")
io.open(MEM, "w", encoding="utf-8").write(m)
n1 = len(m.encode("utf-8"))
print(f"MEMORY.md {n0} -> {n1} B  {'✅ ≤10240' if n1 <= 10240 else f'❌ 超 {n1-10240}'}")
