"""看单 penta6 探针的完整错误。"""
import subprocess
from pathlib import Path

TMP = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\penta_probe")
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
d = TMP / "A"
r = subprocess.run([str(FEBIO), "-i", "t.feb"], cwd=d,
                   capture_output=True, text=True, timeout=120)
out = (r.stdout or "") + (r.stderr or "")
print("returncode", r.returncode)
print("=" * 70)
print(out[-2500:])
print("=" * 70)
lg = d / "t.log"
if lg.exists():
    print("--- t.log ---")
    print(lg.read_text(encoding="utf-8", errors="ignore")[-1500:])
