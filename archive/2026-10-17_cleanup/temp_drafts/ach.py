import io, json

d = json.load(io.open(r'results/opensim_fe/s3_achilles_sweep.json', encoding='utf-8'))
print(type(d))
if isinstance(d, dict):
    for k, v in d.items():
        if isinstance(v, list):
            print(f'{k}: list[{len(v)}]')
        else:
            print(f'{k}: {v}')
    rows = None
    for k, v in d.items():
        if isinstance(v, list) and v and isinstance(v[0], dict):
            rows = v
            print('\n--- rows in', k, '---')
            keys = list(rows[0].keys())
            print('  keys:', keys)
            for r in rows:
                print('  ', {kk: r[kk] for kk in keys[:9]})
            break
elif isinstance(d, list):
    print('rows', len(d), 'keys', list(d[0].keys()))
    for r in d:
        print('  ', r)
