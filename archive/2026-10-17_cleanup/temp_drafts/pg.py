import io

L = io.open(r'docs/reports/总体进度.md', encoding='utf-8').read().split('\n')
for i, l in enumerate(L, 1):
    if ('罚' in l or '1059' in l or '踝子线' in l or 'tied' in l
            or 'k=49' in l or '880' in l):
        print('%4d| %s' % (i, l[:210]))
