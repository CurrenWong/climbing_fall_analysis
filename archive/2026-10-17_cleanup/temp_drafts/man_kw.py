import io

P = r'D:\Project\climbing_fall_analysis\temp\thums\febio_manual.txt'
t = io.open(P, encoding='utf-8', errors='replace').read()
L = t.split('\n')
print(f'manual bytes={len(t.encode("utf-8")):,}  lines={len(L):,}\n')

KW = ['underconstrain', 'under-constrain', 'rigid body mode', 'rigid-body',
      'rigid body', 'initial separation', 'stabiliz', 'stabilis',
      'symmetric_stiffness', 'augmented lagrangian', 'augmented_lagrangian',
      'search_tol', 'node_reloc', 'auto_penalty', 'two_pass',
      'prescribed displacement', 'min_residual', 'velocity', 'quasi-newton',
      'max_ups', 'eigenvalue', 'buckling', 'relax', 'line search']
for k in KW:
    hits = [i + 1 for i, l in enumerate(L) if k in l.lower()]
    print(f'  {k:26s} {len(hits):5d}  ' +
          (f'e.g. L{hits[0]}' if hits else ''))
