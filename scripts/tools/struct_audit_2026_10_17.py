"""结构审计（只读）：① 代码是否按功能分目录 ② 潜在易犯错误（footgun）扫描。

产物：`results/STRUCT_AUDIT_2026-10-17.md`
"""
import collections
import io
import os
import re
import time

BASE = r'D:\Project\climbing_fall_analysis'
SKIP = {'.git', '.venv', '__pycache__', '.pytest_cache', '.ruff_cache', 'node_modules', 'archive'}
TEXT = ('.py', '.md', '.json', '.txt', '.csv')
out = []


def rel(p):
    return os.path.relpath(p, BASE).replace('\\', '/')


def add(s=''):
    print(s)
    out.append(s)


# ---------- ① 目录结构与功能归属 ----------
add("=" * 94)
add("① 目录结构（文件数 / 功能归属）")
add("=" * 94)
add(f"\n{'目录':<34}{'文件数':>7}   功能（按内容判读）")
GUESS = {
    'src/climbing': '**核心库**（1D 两大分支 + 耦合层）',
    'src/climbing/coupling': '**FE 耦合层**（OpenSim→FEBio 链、接触、材料）',
    'scripts/opensim_fe': '**FEBio/OpenSim 探针与驱动**（62）',
    'scripts/ankle_fe': '**踝子线专用**（建模/网格/验收）',
    'scripts/p5_cond': '**论文数字校验**（5 个，成套）',
    'scripts': '⚠️ **根目录是"平的"**：1D 落地线探针 + 杂项工具混放',
    'tests': '**回归测试**（按模块一一对应，命名规约清晰）',
    'docs/methods': '预注册/设计（按 G*/P* 编号，清晰）',
    'docs/reports': '结论/交接/进度（含枢纽 `总体进度.md`）',
    'docs/notes': '踩坑手册 + 旧 TODO（手册 3 本，活跃）',
    'docs/paper': '论文正文（v2）',
    'docs/s5-contact': '⚠️ **S5 阶段遗留**（16 个，10-07）',
    'docs/archive': '已归档文档',
    'results': '产物（报告 + json）',
    'results/opensim_fe': '⚠️ **最大杂堆**（FE 产物 300+）',
    'temp': '⚠️ **草稿区**，但 `temp/opensim_fe` 有 **12,593** 个 FE 产物',
    'model': 'FE 模型输入',
    'paper': '论文源（LaTeX/图表）',
    'video': '短视频线（独立）',
    'legacy': '历史遗留',
    'archive': '本次整理归档',
}
for d in ('src/climbing', 'src/climbing/coupling', 'scripts', 'scripts/opensim_fe',
          'scripts/ankle_fe', 'scripts/p5_cond', 'tests', 'docs/methods', 'docs/notes',
          'docs/reports', 'docs/paper', 'docs/s5-contact', 'docs/archive', 'results',
          'results/opensim_fe', 'temp', 'temp/opensim_fe', 'model', 'paper', 'video',
          'legacy', 'archive'):
    p = os.path.join(BASE, d.replace('/', os.sep))
    if not os.path.isdir(p):
        continue
    n = 0
    for dp, dns, fs in os.walk(p):
        dns[:] = [x for x in dns if x not in SKIP]
        n += len(fs)
    add(f"{d:<34}{n:>7}   {GUESS.get(d, '')}")

# ---------- 收集 py 源码 ----------
srcs = {}
for sub in ('src', 'scripts', 'tests'):
    for dp, dns, fs in os.walk(os.path.join(BASE, sub)):
        dns[:] = [x for x in dns if x not in SKIP]
        for f in fs:
            if f.endswith('.py'):
                p = os.path.join(dp, f)
                try:
                    srcs[rel(p)] = io.open(p, encoding='utf-8', errors='ignore').read()
                except OSError:
                    pass

add("\n" + "=" * 94)
add("② 跨文件重复函数名（同一天里「复制公式」的温床）")
add("=" * 94)
fn = collections.defaultdict(list)
for r, t in srcs.items():
    for m in re.finditer(r'^\s*def (\w+)\(', t, re.M):
        fn[m.group(1)].append(r)
dup = {k: v for k, v in fn.items() if len(set(v)) >= 4 and not k.startswith('_')}
for k, v in sorted(dup.items(), key=lambda kv: -len(set(kv[1])))[:14]:
    add(f"  {k:<32} {len(set(v)):>3} 文件  例: {', '.join(sorted(set(v))[:3])}")

add("\n" + "=" * 94)
add("③ ★ 全局可变「门」（gate dict）—— 隐式状态，脚本间会互相影响")
add("=" * 94)
gates = []
for r, t in srcs.items():
    if not r.startswith('src/climbing/'):
        continue
    for m in re.finditer(r'^([A-Z][A-Z0-9_]*)\s*:\s*dict\s*=\s*\{', t, re.M):
        gates.append((r, m.group(1), t.count(m.group(1))))
for g in gates:
    add(f"  {g[1]:<26} 定义于 {g[0]}")
add("\n  谁在改它们（脚本侧 write 次数）：")
for r, t in srcs.items():
    if not r.startswith('scripts/'):
        continue
    hits = sorted({g[1] for g in gates for _ in [0] if g[1] in t})
    if hits:
        add(f"    {r:<52} {', '.join(hits)}")

add("\n" + "=" * 94)
add("④ 已知坑类的复发检测（模式扫描）")
add("=" * 94)
PATS = [
    (r'np\.gradient\(', '对采样轨迹做数值微分（采样依赖）'),
    (r'parents\[2\]', '三层目录脚本（与 tests 的 parents[2] 不同，易混）'),
    (r'sys\.path\.insert', '脚本自注入 sys.path（须与本文件层级匹配）'),
    (r'_mt\[:, 0\]|m_theta_table', 'M(θ) 表接口（契约见 CANONICAL.md §0）'),
    (r'f_flex_peak_n', '腿力观测量（须与 rhs 同源）'),
    (r'Posture\(|POSTURES\[', '按「键」判定门的姿势（传对象会静默绕门）'),
    (r'to_hdf5|xplt_shim', 'pyfebio 接触 .xplt 无名域缺陷'),
]
for pat, note in PATS:
    hits = [(r, i + 1) for r, t in srcs.items()
            for i, l in enumerate(t.splitlines()) if re.search(pat, l)]
    add(f"\n  ▸ {note}：**{len(hits)} 处**")
    for r, n in sorted(hits)[:6]:
        add(f"      {r}:{n}")
    if len(hits) > 6:
        add(f"      …另 {len(hits) - 6} 处")

io.open(os.path.join(BASE, 'results', 'STRUCT_AUDIT_2026-10-17.md'), 'w',
        encoding='utf-8').write('\n'.join(out) + '\n')
print("\nOUT = results/STRUCT_AUDIT_2026-10-17.md")
