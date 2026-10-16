import io, re

P = r'D:\Project\climbing_fall_analysis\temp\thums\febio_manual.txt'
L = io.open(P, encoding='utf-8', errors='replace').read().split('\n')

TARGETS = [
    'Rigid body modes', 'Forcing convergence', 'Initial Separation',
    'Force-driven Problems', 'Failure to converge', 'Problems due to Contact',
    'No loads applied', 'The penalty method', 'Augmented Lagrangian Method',
    'False Transient Analysis', 'Constraint enforcement', 'Mesh convergence',
    'Geometrical instabilities',
]
for k in TARGETS:
    hits = [i + 1 for i, l in enumerate(L)
            if k.lower() in l.lower() and len(l.strip()) < 90 and i > 1000]
    print(f'{k:28s} -> {hits[:8]}')
