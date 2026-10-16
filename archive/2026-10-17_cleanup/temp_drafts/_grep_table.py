import io
import os

BASE = r'D:\Project\climbing_fall_analysis'
out = []
for d in ('scripts', 'src', 'tests', 'temp'):
    for dp, _, fs in os.walk(os.path.join(BASE, d)):
        for f in fs:
            if not f.endswith('.py'):
                continue
            p = os.path.join(dp, f)
            try:
                txt = io.open(p, encoding='utf-8', errors='ignore').read()
            except OSError:
                continue
            for i, l in enumerate(txt.splitlines()):
                if 'm_theta_table' in l:
                    out.append((os.path.relpath(p, BASE), i + 1, l.strip()))
for t in sorted(out):
    print('%s:%d: %s' % t)
