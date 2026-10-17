"""找 Mesh 段里 discrete_set 的定义语法（2 节点弹簧单元的 elem 写法）。"""
import re

import pdfplumber

PDF = r"D:\Program\FEBioStudio\doc\FEBio_User_Manual.pdf"
with pdfplumber.open(PDF) as pdf:
    pages = [(i + 1, pg.extract_text() or "") for i, pg in enumerate(pdf.pages)]
    for key in ("DiscreteElement", "discrete element set", "DiscreteSet",
                "discrete_set", 'type="discrete"', "<discrete", "set=\"discrete\""):
        hit = [p for p, t in pages if key in t]
        print(f"  {key!r:26} -> 页 {hit[:16]}")
    print("\n--- 含 'Discrete' 且含 '<elem' 的页 ---")
    for p, t in pages:
        if "Discrete" in t and "<elem" in t:
            print(f"  第 {p} 页")
