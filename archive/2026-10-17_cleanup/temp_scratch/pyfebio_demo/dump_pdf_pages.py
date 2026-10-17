"""定向抽取 PDF 指定页的文本。用法: python dump_pdf_pages.py 418 419 420 ..."""
import sys

import pdfplumber

PDF = r"D:\Program\FEBioStudio\doc\FEBio_User_Manual.pdf"
pages = [int(a) for a in sys.argv[1:]] or [418, 419, 420]
with pdfplumber.open(PDF) as pdf:
    n = len(pdf.pages)
    print(f"[PDF 共 {n} 页]")
    for p in pages:
        if p < 1 or p > n:
            print(f"\n### 第 {p} 页 超出范围")
            continue
        t = pdf.pages[p - 1].extract_text() or ""
        print(f"\n{'='*74}\n### 第 {p} 页\n{'='*74}")
        print(t[:4200])
