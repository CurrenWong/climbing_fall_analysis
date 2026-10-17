"""精确定位残留行 + MEMORY 尾部实际文本。"""
import io

R = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"
for f, nm, pats in (
        (R + r"\ANKLE_FIX_TODO.md", "ANKLE_FIX_TODO", ("实际是外翻", "是反的")),
        (R + r"\docs\踝FE可用模型使用说明.md", "使用说明", ("−IE 才是内翻",)),
):
    t = io.open(f, encoding="utf-8").read()
    print(f"===== {nm}")
    for i, l in enumerate(t.splitlines()):
        if any(p in l for p in pats):
            print(f"  {i+1:>4}: {l.strip()[:150]}")

print("\n===== MEMORY 尾部 30 行")
m = io.open(W + r"\MEMORY.md", encoding="utf-8").read()
ls = m.splitlines()
for i, l in enumerate(ls[-30:], start=len(ls) - 30 + 1):
    print(f"  {i:>4}: {l}")
