import io, re, glob, os

W = r'D:\Project\climbing_fall_analysis\temp\opensim_fe\v1_combo'
for d in sorted(glob.glob(os.path.join(W, '*'))):
    lg = glob.glob(os.path.join(d, '*.log'))
    if not lg:
        continue
    txt = io.open(lg[0], encoding='utf-8', errors='replace').read()
    neg = [int(m) for m in re.findall(r'(\d+)\s+negative jacobians detected', txt)]
    times = [float(m) for m in re.findall(r'failed to converge at time\s*:\s*([\d.eE+-]+)', txt)]
    steps = re.findall(r'beginning time step \d+ : ([\d.eE+-]+)', txt)
    print(f'=== {os.path.basename(d)} ===')
    print(f'  负 Jacobian 次数 {len(neg)}  计数: min={min(neg) if neg else "-"} '
          f'max={max(neg) if neg else "-"} last={neg[-3:] if neg else "-"}')
    print(f'  失败时刻(最后 5): {times[-5:]}')
    print(f'  进入的时间步数: {len(steps)}  最后进入 t={steps[-1] if steps else "-"}')
    print()
