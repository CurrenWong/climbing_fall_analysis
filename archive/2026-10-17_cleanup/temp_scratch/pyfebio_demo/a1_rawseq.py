from pathlib import Path

W = Path(r"D:\Project\climbing_fall_analysis\temp/pyfebio_demo/a1_matrix")
for tag in ("talus", "calcaneus"):
    s = (W / f"{tag}.febio.log").read_text(encoding="utf-8", errors="replace")
    i = s.find("beginning time step 1")
    seg = [l for l in s[i:i + 3000].split("\n")][:52]
    print("=" * 10, tag)
    print("\n".join(seg))
