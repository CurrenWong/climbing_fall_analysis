"""找「rigid body」材料的定义语法。"""
import re

import pdfplumber

PDF = r"D:\Program\FEBioStudio\doc\FEBio_User_Manual.pdf"
with pdfplumber.open(PDF) as pdf:
    pages = [(i + 1, pg.extract_text() or "") for i, pg in enumerate(pdf.pages)]
for p, t in pages:
    for m in re.finditer(r'<material[^>]*type\s*=\s*"(rigid[^"]*)"', t):
        print(f"  第 {p} 页: {m.group(0)}")
print("\n--- 含 'rigid body material' 或 'Rigid Material' 的页 ---")
for p, t in pages:
    if re.search(r"rigid (body )?material", t, re.I):
        print(f"  第 {p} 页", "|", " ".join(
            l.strip() for l in t.splitlines()
            if re.search(r"rigid (body )?material", l, re.I))[:150])
