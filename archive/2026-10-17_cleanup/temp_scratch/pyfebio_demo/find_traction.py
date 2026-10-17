"""定位 surface_load / traction 的 lc 用法。"""
import re

import pdfplumber

PDF = r"D:\Program\FEBioStudio\doc\FEBio_User_Manual.pdf"
with pdfplumber.open(PDF) as pdf:
    pages = [(i + 1, pg.extract_text() or "") for i, pg in enumerate(pdf.pages)]
for k in ("surface_load", "Surface Loads", "Surface Traction", "traction"):
    hit = [p for p, t in pages if k in t]
    print(f"{k!r:20} -> 页 {hit[:16]}")
print("\n--- 含 surface_load 的页上，traction 附近文本 ---")
for p, t in pages:
    if "surface_load" in t and "traction" in t:
        for m in re.finditer(r"traction", t):
            a, b = max(0, m.start() - 500), min(len(t), m.end() + 500)
            print(f"\n### p{p}")
            print(t[a:b].replace("\n", " | "))
            break
        break
