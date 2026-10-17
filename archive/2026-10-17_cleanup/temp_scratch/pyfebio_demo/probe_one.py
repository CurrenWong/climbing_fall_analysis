import subprocess
from pathlib import Path

d = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\penta_probe\A")
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
r = subprocess.run([str(FEBIO), "-i", "t.feb"], cwd=d,
                   capture_output=True, text=True, timeout=180)
print("rc", r.returncode)
print("--- stdout tail ---")
print((r.stdout or "")[-1800:])
print("--- stderr ---")
print((r.stderr or "")[-800:])
lg = d / "t.log"
print("--- log exists:", lg.exists())
if lg.exists():
    print(lg.read_text(encoding="utf-8", errors="ignore")[-1200:])
