from pathlib import Path

t = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\_pm1run.txt").read_text(
    encoding="utf-8", errors="replace")
i = t.find("Negative jacobian detected.")
print("--- 报错前 1800 字符 ---")
print(t[max(0, i - 1800):i + 300])
