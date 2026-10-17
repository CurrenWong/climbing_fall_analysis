"""核对生成的 .feb：最大节点号 / ground_rigid_ns / rigid BC 是否引用了不存在的节点。"""
import re
from pathlib import Path

F = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio_lig\ankle_lig_WK.feb")
t = F.read_text(encoding="latin-1", errors="replace")
ids = [int(m.group(1)) for m in re.finditer(r'<node id="(\d+)">', t)]
print(f"节点数 {len(ids)}  最大 id {max(ids)}  是否连续 1..max: {ids == list(range(1, max(ids)+1))}")
for nm in ("talus_rigid_ns", "ground_rigid_ns", "talus_low_ns"):
    m = re.search(rf'<NodeSet name="{nm}">([^<]*)</NodeSet>', t)
    if not m:
        print(f"  {nm}: (无)")
        continue
    v = [int(x) for x in m.group(1).split(",") if x.strip()]
    bad = [x for x in v if x < 1 or x > max(ids)]
    print(f"  {nm}: {len(v)} 个节点，越界 {len(bad)} 个 {bad[:8]}")
print("\n--- Rigid 段 ---")
i = t.find("<Rigid>")
j = t.find("</Rigid>")
print(t[i:j + 9] if i >= 0 else "(无 <Rigid>)")
print("\n--- rigid BC ---")
for m in re.finditer(r'<bc name="(rigid_[^"]+)"[^>]*>(.*?)</bc>', t, re.S):
    print(f"  {m.group(1)}: {m.group(2).strip()[:120]}")
