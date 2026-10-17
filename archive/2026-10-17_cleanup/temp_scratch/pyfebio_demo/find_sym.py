"""查 symmetric_stiffness 在 FEBio 4.13 手册里的正确位置。"""
import re

import pdfplumber

PDF = r"D:\Program\FEBioStudio\doc\FEBio_User_Manual.pdf"
with pdfplumber.open(PDF) as pdf:
    pages = [(i + 1, pg.extract_text() or "") for i, pg in enumerate(pdf.pages)]
for k in ("symmetric_stiffness", "symmetric stiffness"):
    hit = [p for p, t in pages if k in t]
    print(f"{k!r} -> 页 {hit[:14]}")
for p, t in pages:
    if "symmetric_stiffness" in t:
        for m in re.finditer(r"symmetric_stiffness", t):
            a, b = max(0, m.start() - 700), min(len(t), m.end() + 400)
            print(f"\n--- p{p} ---")
            print(t[a:b].replace("\n", " | "))
        break
