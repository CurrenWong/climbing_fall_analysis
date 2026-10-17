"""在 FEBio User Manual 里定位 Rigid 段语法（刚体 + reference node + 转动约束）。"""
import re

import pdfplumber

PDF = r"D:\Program\FEBioStudio\doc\FEBio_User_Manual.pdf"
KEYS = ("rigid_body", "RigidDomain", "rigid_constraint", "rigid_bc", "rigid_load",
        "reference node", "rigid material", 'type="rigid"', "Rigid Section",
        "rigid body", "RigidBody")
with pdfplumber.open(PDF) as pdf:
    pages = [(i + 1, pg.extract_text() or "") for i, pg in enumerate(pdf.pages)]
print(f"[PDF {len(pages)} 页]")
for k in KEYS:
    hit = [p for p, t in pages if k.lower() in t.lower()]
    print(f"  {k!r:22} -> 页 {hit[:16]}")
print("\n--- 含 'Rigid' 的章节标题页 ---")
for p, t in pages:
    for m in re.finditer(r"^\s*(\d+(?:\.\d+)*)\s+Rigid.*$", t, re.M):
        print(f"  第 {p} 页: {m.group(0).strip()[:80]}")
