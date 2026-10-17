"""从 FEBio User Manual 抽取 tied 相关章节原文（pdfplumber）。"""
import re

import pdfplumber

PDF = r"D:\Program\FEBioStudio\doc\FEBio_User_Manual.pdf"
PAT = re.compile(r"tied|TiedInterface|tie contact", re.I)

hits = []
with pdfplumber.open(PDF) as pdf:
    for i, pg in enumerate(pdf.pages):
        t = pg.extract_text() or ""
        if PAT.search(t):
            hits.append((i + 1, t))

print(f"总页 {i + 1} / 命中 {len(hits)} 页")
print("=" * 74)
for n, t in hits:
    lines = t.split("\n")
    for j, ln in enumerate(lines):
        if PAT.search(ln):
            lo, hi = max(0, j - 4), min(len(lines), j + 10)
            print(f"\n--- p.{n} ---")
            print("\n".join(lines[lo:hi]))
            break
