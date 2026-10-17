"""定位 FEBio 手册里的 Step/多步分析语法章节。"""
import re

import pdfplumber

PDF = r"D:\Program\FEBioStudio\doc\FEBio_User_Manual.pdf"
with pdfplumber.open(PDF) as pdf:
    pages = [(i + 1, pg.extract_text() or "") for i, pg in enumerate(pdf.pages)]
print(f"[PDF {len(pages)} 页]")
print("--- 章节标题含 Step / Multi-step / Analysis ---")
for p, t in pages:
    for m in re.finditer(r"^\s*(\d+(?:\.\d+)*)\s+(Step|Multi-?step|Analysis|"
                         r"Boundary|Loads)\b.*$", t, re.M):
        s = m.group(0).strip()
        if len(s) < 70:
            print(f"  p{p:>4}: {s}")
    if p > 230:
        break
print("\n--- 关键词命中 ---")
for k in ("<Step", "</Step>", "<step id=", "multi-step", "Multi-step analysis",
          "module=\"solid\"", "<Boundary>", "<Loads>"):
    hit = [p for p, t in pages if k in t]
    print(f"  {k!r:22} -> 页 {hit[:14]}")
