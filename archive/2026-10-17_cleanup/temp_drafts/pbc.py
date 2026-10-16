import io, json

R = r'D:\Project\climbing_fall_analysis\results\opensim_fe'
d = json.load(io.open(R + r'\plantar_bc_summary.json', encoding='utf-8'))
print('meta:', json.dumps(d['meta'], ensure_ascii=False)[:600])
print('load:', json.dumps(d['load'], ensure_ascii=False)[:300])
print('\ncases:')
for c in d['cases']:
    keep = {k: c[k] for k in list(c.keys())[:12]}
    print('  ', json.dumps(keep, ensure_ascii=False)[:300])
print('\nspring_scan:')
ss = d['spring_scan']
print('  type', type(ss), len(ss) if hasattr(ss, '__len__') else '')
if isinstance(ss, list):
    for r in ss[:10]:
        print('  ', json.dumps(r, ensure_ascii=False)[:240])
elif isinstance(ss, dict):
    print('  keys', list(ss.keys())[:10])
    print('  ', json.dumps(ss, ensure_ascii=False)[:800])
print('\nfirst_fracture_extrapolation:',
      json.dumps(d['first_fracture_extrapolation'], ensure_ascii=False)[:600])
