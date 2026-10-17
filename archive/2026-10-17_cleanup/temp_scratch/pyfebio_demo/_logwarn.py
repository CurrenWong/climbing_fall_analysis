import io
from pathlib import Path

t = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio"
         r"\ankle_tibiotalar.log").read_text(encoding="latin-1", errors="replace")
L = t.split("\n")
print("\n".join(f"{i+1}: {l.rstrip()}" for i, l in enumerate(L[:80]) if l.strip()))
