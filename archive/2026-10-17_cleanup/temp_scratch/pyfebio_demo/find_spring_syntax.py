"""从 FEBio User Manual PDF 里找「弹簧边界条件」的语法。

优先 pdfplumber；没有就退回 pypdf/PyPDF2。
搜索关键词：type="spring"、spring_stiffness、<bc ... spring
"""
import re
import sys

PDF = r"D:\Program\FEBioStudio\doc\FEBio_User_Manual.pdf"
KEYS = ('type="spring"', "spring_stiffness", "spring_", "linear spring",
        "elastic spring", "prescribed displacement")


def pages_text():
    try:
        import pdfplumber
        print("[extractor] pdfplumber")
        with pdfplumber.open(PDF) as pdf:
            for i, pg in enumerate(pdf.pages):
                yield i + 1, (pg.extract_text() or "")
        return
    except ImportError:
        pass
    try:
        from pypdf import PdfReader
        print("[extractor] pypdf")
    except ImportError:
        from PyPDF2 import PdfReader
        print("[extractor] PyPDF2")
    for i, pg in enumerate(PdfReader(PDF).pages):
        yield i + 1, (pg.extract_text() or "")


hits = []
for pno, txt in pages_text():
    for k in KEYS:
        if k in txt:
            hits.append((pno, k))
            break

print(f"命中页数 {len(hits)}: {[h[0] for h in hits][:40]}")
# 打印最可能的几页（含 type="spring" 或 spring_stiffness）
want = [p for p, k in hits if k in ('type="spring"', "spring_stiffness")]
print(f"强命中的页: {want[:20]}")
for pno, txt in pages_text():
    if pno in want[:3]:
        print(f"\n{'='*74}\n### PDF 第 {pno} 页\n{'='*74}")
        print(txt[:3000])
