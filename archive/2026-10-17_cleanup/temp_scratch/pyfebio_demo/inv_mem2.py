from pathlib import Path

WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
p = WS / "memory/2026-10-04.md"
print("=" * 70)
print(p.name, len(p.read_text(encoding="utf-8").encode("utf-8")), "bytes")
print(p.read_text(encoding="utf-8")[:700])
for f in ("memory/2026-10-04/ankle-fe-solve-and-handover.md",):
    q = WS / f
    if q.exists():
        s = q.read_text(encoding="utf-8")
        print("=" * 70)
        print(f, len(s.encode("utf-8")), "bytes")
        for l in s.split("\n"):
            if l.startswith("#"):
                print("   ", l.strip()[:66])
