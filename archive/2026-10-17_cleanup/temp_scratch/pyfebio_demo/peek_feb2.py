"""看 .feb 的元素段结构。"""
import re
from pathlib import Path

F = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio_lig\ankle_lig_WK.feb")
t = F.read_text(encoding="utf-8")
print("总长", len(t))
# 所有一级标签
print("\n顶层标签:")
for m in re.finditer(r"^  <([A-Za-z_]+)", t, re.M):
    print("   ", m.group(1))
i = t.find("<Mesh>")
print("\n--- <Mesh> 开头 700 字符 ---")
print(t[i:i + 700] if i >= 0 else "无 <Mesh>")
j = t.find("tet4")
print("\n--- 第一处 tet4 附近 500 字符 ---")
print(t[max(0, j - 200):j + 300] if j >= 0 else "无 tet4")
