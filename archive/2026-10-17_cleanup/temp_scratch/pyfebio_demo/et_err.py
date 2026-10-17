import subprocess
from pathlib import Path

d = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\elem_test\penta6")
r = subprocess.run([r"D:\Program\FEBioStudio\bin\febio4.exe", "-i", "t.feb"],
                   cwd=d, capture_output=True, text=True, timeout=180)
o = (r.stdout or "") + (r.stderr or "")
i = o.find("FAILED")
print(o[max(0, i - 300):i + 600] if i > 0 else o[-900:])
print("=== files ===")
for f in sorted(d.iterdir()):
    print("  ", f.name)
