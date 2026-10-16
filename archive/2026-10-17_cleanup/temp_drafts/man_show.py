import io
import sys

P = r'D:\Project\climbing_fall_analysis\temp\thums\febio_manual.txt'
t = io.open(P, encoding='utf-8', errors='replace').read()
L = t.split('\n')


def show(a, b, tag=''):
    print(f'===== lines {a}-{b} {tag} =====')
    for i in range(a - 1, min(b, len(L))):
        print(f'{i+1:6d}| {L[i][:190]}')
    print()


if __name__ == '__main__':
    a, b = int(sys.argv[1]), int(sys.argv[2])
    show(a, b, sys.argv[3] if len(sys.argv) > 3 else '')
