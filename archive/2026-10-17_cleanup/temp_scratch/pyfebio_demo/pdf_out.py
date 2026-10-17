"""查手册 Output 段：plotfile 类型 + logfile/node_data（拿文本读数的最省事路径）。"""
import re

import pdfplumber

PDF = r"D:\Program\FEBioStudio\doc\FEBio_User_Manual.pdf"


def dump(a, b, title):
    print(f"\n{'='*76}\n=== {title}  p{a}~{b} ===\n{'='*76}")
    with pdfplumber.open(PDF) as pdf:
        for i in range(a, min(b + 1, len(pdf.pages) + 1)):
            t = pdf.pages[i - 1].extract_text() or ""
            print(f"\n----- p{i} -----")
            print(t[:2600])


# 找 Output Section 的页码
with pdfplumber.open(PDF) as pdf:
    for i, pg in enumerate(pdf.pages[:14], 1):
        t = pg.extract_text() or ""
        for l in t.split("\n"):
            if re.search(r"Output Section|Plotfile|Logfile|Output Variables|Nodal Output",
                         l, re.I):
                print(f"[TOC p{i}] {l.strip()[:110]}")
