import importlib.util as u
for m in ("fitz", "pdfplumber", "PyPDF2", "pypdf", "pymupdf"):
    print(f"{m:12s} {bool(u.find_spec(m))}")
