import subprocess
from pathlib import Path

for name in ("penta6", "tet4x3", "penta15"):
    d = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\elem_test2") / name
    r = subprocess.run([r"D:\Program\FEBioStudio\bin\febio4.exe", "-i", "t.feb"],
                       cwd=d, capture_output=True, text=True, timeout=180)
    o = (r.stdout or "") + (r.stderr or "")
    i = o.find("FAILED")
    print(f"===== {name} =====")
    print(o[max(0, i - 100):i + 400] if i > 0 else o[-600:])
