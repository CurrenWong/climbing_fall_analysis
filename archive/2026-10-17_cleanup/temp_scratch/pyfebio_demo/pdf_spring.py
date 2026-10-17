"""从 FEBio 用户手册里提取 spring 元素的权威 XML 语法。"""
import re

import pdfplumber

PDF = r"D:\Program\FEBioStudio\doc\FEBio_User_Manual.pdf"

pats = [r"linear spring", r"spring_constant", r"<spring", r"Springs"]
with pdfplumber.open(PDF) as pdf:
    print(f"总页数 {len(pdf.pages)}")
    hits = []
    for i, pg in enumerate(pdf.pages, 1):
        try:
            t = pg.extract_text() or ""
        except Exception as e:                                    # noqa: BLE001
            print(f"  page {i} 提取失败 {e}")
            continue
        if re.search(r"linear spring|spring_constant", t, re.I):
            hits.append((i, t))
    print(f"命中 {len(hits)} 页: {[h[0] for h in hits]}\n")

    # 目录里找 Spring 章节
    with pdfplumber.open(PDF) as pdf2:
        for i, pg in enumerate(pdf2.pages, 1):
            t = pg.extract_text() or ""
            if re.search(r"spring", t, re.I) and re.search(r"\.{2,}\s*\d+", t):
                for l in t.split("\n"):
                    if re.search(r"spring", l, re.I):
                        print(f"  [目录 p{i}] {l.strip()[:100]}")
        if len(pdf2.pages) > 12:
            print("\n(目录扫描仅前 12 页)")

    for i, t in hits[:4]:
        print(f"\n{'='*74}\n===== 手册 p{i} =====\n{'='*74}")
        print(t[:3000])
