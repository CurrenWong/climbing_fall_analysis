"""看 .feb 的实际标签格式。"""
import re
from pathlib import Path

F = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio_lig\ankle_lig_WK.feb")
t = F.read_text(encoding="utf-8")
print("总长", len(t))
for tag in ("Nodes", "node", "Surface", "tri3", "quad4", "MeshDomains",
            "Solid", "Material"):
    n = len(re.findall(rf"<{tag}[ >]", t))
    print(f"  <{tag}> {n} 次")

i = t.find("<Nodes>")
print("\n--- Nodes 开头 300 字符 ---")
print(t[i:i + 300])

j = t.find("<Surface")
print("\n--- 第一个 Surface 开头 300 字符 ---")
print(t[j:j + 300])

k = t.find("<MeshDomains>")
print("\n--- MeshDomains 开头 400 字符 ---")
print(t[k:k + 400])

print("\n--- 所有 Surface name ---")
for m in re.finditer(r'<Surface name="([^"]+)"', t):
    print("   ", m.group(1))
