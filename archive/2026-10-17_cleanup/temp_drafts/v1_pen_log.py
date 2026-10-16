import io, glob, os, re

W = r'D:\Project\climbing_fall_analysis\temp\opensim_fe\v1_pen'
for d in sorted(glob.glob(os.path.join(W, '*'))):
    lg = glob.glob(os.path.join(d, '*.log'))
    tag = os.path.basename(d)
    if not lg:
        print(f'{tag:16s} no log')
        continue
    t = io.open(lg[0], encoding='utf-8', errors='replace').read()
    neg = [int(m) for m in re.findall(r'(\d+)\s+negative jacobians detected', t)]
    ts = re.findall(r'failed to converge at time\s*:\s*([\d.eE+-]+)', t)
    unk = re.findall(r'tag "([^"]+)" : unrecognized', t)
    steps = re.findall(r'beginning time step \d+ : ([\d.eE+-]+)', t)
    pen = re.findall(r'<penalty>([^<]+)</penalty>', io.open(
        glob.glob(os.path.join(d, '*.feb'))[0], encoding='utf-8',
        errors='replace').read()) if glob.glob(os.path.join(d, '*.feb')) else []
    print(f'{tag:16s} penalty={pen[0] if pen else "?":>8s}  negJ={len(neg):3d} '
          f'first={neg[0] if neg else "-"} max={max(neg) if neg else "-"}  '
          f'进入步={len(steps)}  最后 t={steps[-1] if steps else "-"}  '
          f'unrecognized={set(unk) if unk else "-"}')
