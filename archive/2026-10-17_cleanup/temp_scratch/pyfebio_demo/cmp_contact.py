"""对比两个 .feb 的 <Contact> 段 + Surface 定义（找出差异）。"""
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
FILES = {
    "lig (我, 新格式+TIE)": ROOT / "temp/thums/febio_lig/ankle_lig_WK.feb",
}
# 找 build_feb.py 的输出
for p in sorted((ROOT / "temp/thums").rglob("*.feb")):
    if "lig" not in p.name.lower():
        FILES[f"build_feb: {p.name}"] = p


def section(txt, tag):
    i = txt.find(f"<{tag}>")
    if i < 0:
        i = txt.find(f"<{tag} ")
    if i < 0:
        return f"(无 <{tag}>)"
    j = txt.find(f"</{tag}>", i)
    return txt[i:j + len(tag) + 3] if j > 0 else txt[i:i + 2500]


for name, p in FILES.items():
    if not p.exists():
        print(f"(缺) {name}")
        continue
    t = p.read_text(encoding="utf-8", errors="replace")
    print(f"\n{'='*78}\n### {name}   ({p.stat().st_size/1024:.0f} KB)\n{'='*78}")
    c = section(t, "Contact")
    print(f"--- Contact ({len(c)} chars) ---")
    print(c[:2200])
