"""查 FEBio 对 lc 的警告 + dump LoadData/traction 实际内容。"""
import io
import re

t = io.open(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\_ms_febio.txt",
            encoding="utf-8", errors="replace").read()
print("--- 含 lc/load/traction/attribute/invalid 的行 ---")
for l in t.splitlines():
    ls = l.strip()
    if ls and any(k in ls.lower() for k in
                  ("lc=", "loadcurve", "load_controller", "traction", "ignor",
                   "not supported", "unknown", "attribute", "invalid", "warn")):
        print("  ", ls[:100])

f = io.open(r"D:\Project\climbing_fall_analysis\temp\thums\febio_lig"
            r"\ankle_lig_WK.feb", encoding="latin-1").read()
i, j = f.find("<LoadData>"), f.find("</LoadData>")
print("\n===== <LoadData> =====")
print(f[i:j + 11])
k = f.find("<surface_load")
print("===== step2 的 <Loads> =====")
print(f[k:k + 300])
