"""修 daily 索引：守卫用 LRN ID（不再用节号，避免 §十七 撞号静默跳过）。"""
from pathlib import Path

D = Path(r"C:\Users\curre\.copaw\workspaces\default\memory\2026-10-04.md")
L = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\_index_line.md")
MARK = "LRN-20261005-059"

t = D.read_text(encoding="utf-8")
line = L.read_text(encoding="utf-8").strip()
if MARK not in t:
    D.write_text(t.rstrip() + "\n" + line + "\n", encoding="utf-8")
    print(f"OK daily 已追加 -> {D.stat().st_size:,} B")
else:
    print("daily 已含该索引，跳过")
tt = D.read_text(encoding="utf-8")
print(f"  [{'OK' if MARK in tt else 'XX'}] daily 含 {MARK}")
print(f"  [{'OK' if 'dispx' in tt else 'XX'}] daily 含 dispx")
print(f"  [{'OK' if 'acceptance.py' in tt else 'XX'}] daily 含 acceptance.py")
