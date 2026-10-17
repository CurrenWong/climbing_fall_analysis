"""找附录 E 里 node_data / model_data 的合法变量名（要 ux/uy/uz、反力、|u| 等）。"""
import re

import pdfplumber

PDF = r"D:\Program\FEBioStudio\doc\FEBio_User_Manual.pdf"

with pdfplumber.open(PDF) as pdf:
    n = len(pdf.pages)
    # 附录 E 通常在末尾 40 页内
    for i in range(max(1, n - 45), n + 1):
        t = pdf.pages[i - 1].extract_text() or ""
        if re.search(r"\bux\b.*\buy\b|\bRx\b.*\bRy\b|Nodal data|Data expression|log variable",
                     t, re.I):
            print(f"\n{'='*76}\n=== p{i} ===\n{'='*76}")
            print(t[:3000])
