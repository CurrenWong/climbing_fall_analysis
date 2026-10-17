"""从 Anderson 2006 原文 PDF 里抓 §2.2 边界条件相关段落（spring / rigid / restraint）。"""
import re

import pdfplumber

PDF = r"D:\Project\climbing_fall_analysis\paper\2006FAIankleligaments.pdf"
PAT = re.compile(r"(spring|rigid|restrain|restraint|inversion|internal-external|"
                 r"seating|preload|pre-load)", re.I)

with pdfplumber.open(PDF) as pdf:
    pages = [(i + 1, pg.extract_text() or "") for i, pg in enumerate(pdf.pages)]
print(f"[PDF {len(pages)} 页]")
for p, t in pages:
    if "spring" in t.lower():
        hits = sorted({m.group(0).lower() for m in PAT.finditer(t)})
        print(f"  第 {p} 页  命中 {hits}")

# 打印含 spring 的那几页的 spring 上下文
print("\n" + "=" * 78)
for p, t in pages:
    if "spring" in t.lower():
        print(f"\n########## 第 {p} 页 ##########")
        for m in re.finditer(r"spring", t, re.I):
            a, b = max(0, m.start() - 900), min(len(t), m.end() + 900)
            print("..." + t[a:b].replace("\n", " ") + "...")
            print("-" * 70)
        print()
