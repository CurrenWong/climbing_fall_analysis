"""查手册：Boundary 段的 spring BC 语法 + Mesh 段 discrete_set 语法。"""
import re

import pdfplumber

PDF = r"D:\Program\FEBioStudio\doc\FEBio_User_Manual.pdf"
KEYS = [r"type=\"spring\"", r"spring_constant", r"discrete_set", r"discrete element set",
        r"<Discrete>", r"Linear Spring", r"Spring\b"]

with pdfplumber.open(PDF) as pdf:
    for i, pg in enumerate(pdf.pages, 1):
        t = pg.extract_text() or ""
        if not re.search(r"spring", t, re.I):
            continue
        # 只关心含 XML 片段或 parameter 表的页
        if not re.search(r"<|parameter\s+Description", t):
            continue
        tags = [k for k in KEYS if re.search(k, t, re.I)]
        if not tags:
            continue
        print(f"\n{'='*76}\n=== 手册 p{i}  命中 {tags} ===\n{'='*76}")
        # 打印含 spring 的段落 + 其前后文
        lines = t.split("\n")
        show = set()
        for j, l in enumerate(lines):
            if re.search(r"spring", l, re.I):
                for k in range(max(0, j - 6), min(len(lines), j + 12)):
                    show.add(k)
        for j in sorted(show):
            print(f"  {j:>3}|{lines[j]}")
