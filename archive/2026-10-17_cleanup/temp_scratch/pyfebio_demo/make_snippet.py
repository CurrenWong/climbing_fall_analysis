"""回测：把「修前」的死代码片段写出来，确认扫描器抓得住。"""
import io
from pathlib import Path

SNIP = '''for surf, F, _faces in (("tibia_top", F_t, top_tibia),):
    if not _use_force or not _faces or F <= 0:
        continue
        # traction = 力 / 面积
        a = area_proj(top_tibia)
        tr = F / a if a > 0 else 0.0
        w("    <surface_load ...>")
'''
p = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\_snippet_orig.py")
p.write_text(SNIP, encoding="utf-8")
print(f"写入 {p} ({p.stat().st_size} B)")
print("--- 片段 ---")
print(SNIP)
