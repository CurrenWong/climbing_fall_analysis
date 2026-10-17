"""dump FEBio 手册 p.145-150 全文（tied 系列完整参数）。"""
import pdfplumber

PDF = r"D:\Program\FEBioStudio\doc\FEBio_User_Manual.pdf"
with pdfplumber.open(PDF) as pdf:
    for p in range(145, 150):
        t = pdf.pages[p - 1].extract_text() or ""
        print(f"\n{'=' * 30} p.{p} {'=' * 30}")
        print(t)
