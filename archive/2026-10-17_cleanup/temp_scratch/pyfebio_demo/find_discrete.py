"""搜 PDF 里 Discrete 段的完整语法（discrete_material / discrete_set 等）。"""
import re

import pdfplumber

PDF = r"D:\Program\FEBioStudio\doc\FEBio_User_Manual.pdf"
KEYS = ("discrete_material", "discrete_set", "linear spring", "type=\"discrete\"",
        "Discrete section", "spring", "nonlinear spring")

with pdfplumber.open(PDF) as pdf:
    pages = [(i + 1, pg.extract_text() or "") for i, pg in enumerate(pdf.pages)]
    print(f"[PDF {len(pages)} 页]")
    for key in KEYS[:4]:
        hit = [p for p, t in pages if key in t]
        print(f"  {key!r:24} -> 页 {hit[:14]}")
    print()
    # 找 "Discrete" 章节标题页
    for p, t in pages:
        if re.search(r"^\s*\d+(\.\d+)*\s+Discrete", t, re.M) or \
           re.search(r"\nDiscrete\s*\n.*spring", t):
            print(f"  [章节候选] 第 {p} 页: {t.strip().splitlines()[0][:70]}")
