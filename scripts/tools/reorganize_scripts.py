"""`scripts/` 根下按功能分目录（**修正版**）。

首版三个缺陷（2026-10-17 事故，同一类"静默部分失败"）：
  * ① **没建目标目录** ⇒ `git mv` 全部失败，但 `check=False` 不抛错；
  * ② **无条件打印「已移动 35 个」** —— 统计的是"计划"，不是"成功" ⇒ **脚本在说谎**；
  * ③ ② 步没防"文件不存在"（失败的那些）⇒ 抛异常 ⇒ **③ 从未执行**。

本版：先 `makedirs`；**按实际存在性判定**（src 在⇒移，dst 在⇒视为已移）；每步计数**以磁盘为准**；
② 步 `try/except` 跳过缺失；③ 步只在文件真被移动/已在目标时才重写引用。
"""
import argparse
import io
import os
import re
import subprocess
import sys

BASE = r'D:\Project\climbing_fall_analysis'
SKIP = {'.git', '.venv', '__pycache__', '.pytest_cache', '.ruff_cache', 'node_modules', 'archive'}
TEXT_EXT = ('.py', '.md', '.json', '.txt', '.toml', '.cfg', '.yml', '.yaml', '.ini')

MAP = {
    '1d': ['chain_crush_probe', 'chain_soft_probe', 'chain_tissue_probe',
           'ecc_muscle_probe', 'fp_anchor_probe', 'fp_contact_probe',
           'fp_stroke_probe', 'leg2seg_probe', 'pad_energy_check',
           'passive_path_probe', 'second_contact_probe', 'unload_process_probe',
           'recompute_hstar'],
    'tools': ['archive_cleanup', 'archive_doc_clusters', 'cleanup_parallel',
              'cleanup_round2', 'scan_stale_refs', 'audit_2026_10_17',
              'sweep2_2026_10_17', 'struct_audit_2026_10_17', 'smoke',
              'reorganize_scripts'],
    'paper': ['calibrate_rope', 'chain_vtarget_probe', 'mc_stat_robustness',
              'phase_b2_analysis', 'phase_p1_montecarlo', 'rigid_vtarget_probe',
              'validate_vs_fem', 'verify_b3', 'verify_rope', 'wave_rod_probe'],
    'figures': ['make_ascii_flowchart', 'make_flowchart', 'visualize_b2'],
}


def rel(p):
    return os.path.relpath(p, BASE).replace('\\', '/')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--go', action='store_true')
    a = ap.parse_args()

    pairs = []
    for sub, names in MAP.items():
        for n in names:
            pairs.append((f'scripts/{n}.py', f'scripts/{sub}/{n}.py'))
    print(f"计划 {len(pairs)} 个映射")
    if not a.go:
        for old, new in pairs:
            print(f"    {old} → {new}")
        print("(dry-run)")
        return 0

    # ① 移动（**先建目录**；以磁盘为准计数）
    n_moved = n_already = 0
    for old, new in pairs:
        src = os.path.join(BASE, old.replace('/', os.sep))
        dst = os.path.join(BASE, new.replace('/', os.sep))
        if os.path.isfile(src):
            os.makedirs(os.path.dirname(dst), exist_ok=True)   # ★ 首版缺这句
            r = subprocess.run(['git', 'mv', old, new], cwd=BASE,
                               capture_output=True, text=True)
            if r.returncode != 0 or not os.path.isfile(dst):
                print(f"  ❌ git mv 失败 {old} → {new}: {r.stderr.strip()[:120]}")
            else:
                n_moved += 1
        elif os.path.isfile(dst):
            n_already += 1
        else:
            print(f"  ⚠️ 两边都不存在: {old} / {new}")
    print(f"① 实际移动 **{n_moved}** 个；已在目标 {n_already} 个（计划 {len(pairs)}）")

    # ② parents[N] += 1（只对**真在子目录里**的文件）
    fixed = 0
    for old, new in pairs:
        p = os.path.join(BASE, new.replace('/', os.sep))
        if not os.path.isfile(p):
            continue
        s = io.open(p, encoding='utf-8', errors='ignore').read()
        s2 = re.sub(r'parents\[(\d+)\]', lambda m: f'parents[{int(m.group(1)) + 1}]', s)
        if s2 != s:
            io.open(p, 'w', encoding='utf-8').write(s2)
            fixed += 1
    print(f"② parents[N] +1：{fixed} 个文件")

    # ③ 引用重写（只对**移动已生效**的映射）
    live = [(old, new) for old, new in pairs
            if os.path.isfile(os.path.join(BASE, new.replace('/', os.sep)))]
    rew = {}
    for dp, dns, fs in os.walk(BASE):
        dns[:] = [d for d in dns if d not in SKIP]
        for f in fs:
            if not f.endswith(TEXT_EXT):
                continue
            p = os.path.join(dp, f)
            r = rel(p)
            if r.startswith(('archive/', 'temp/')):
                continue
            try:
                s = io.open(p, encoding='utf-8', errors='ignore').read()
            except OSError:
                continue
            o = s
            c = 0
            for old, new in live:
                if old in s:
                    c += s.count(old)
                    s = s.replace(old, new)
            if s != o:
                io.open(p, 'w', encoding='utf-8').write(s)
                rew[r] = c
    print(f"③ 引用重写：{len(rew)} 个文件，共 {sum(rew.values())} 处")
    for r, n in sorted(rew.items())[:15]:
        print(f"    {r}  ({n})")
    return 0


if __name__ == '__main__':
    sys.exit(main())
