# -*- coding: utf-8 -*-
"""两步法探针（手册 §6.1.1 / §6.1.2）—— 「先用位移压到接触，再换成力」。

手册 §8.5.3 原文给的通用正解：

> Another way around this issue is to **first use a displacement to bring the surfaces in
> contact and then replace the displacement with a force.**

手册 §6.1.2 明确多步分析**正是为此设计**：

> …**This makes it possible, for example, to switch between load and displacement
> controlled boundary conditions in multi-step analyses.**

且 §6.1.1 规定：**放在第一个 Step 之前**的 Boundary/Loads/Contact 是 **global**（各步都生效），
放在 Step 内的只对那一步生效。pyfebio 生成的 deck 恰好把这三块都放在 `<Steps>` 之前 ⇒
**global**。故本探针做最小 XML 手术：

1. **Step 1**（新增，插在最前）：`<Boundary><prescribe bc="y" lc=2>` 把 `subtalar_joint`
   节点沿 −y 压下 `pretouch_mm`（> 初始间隙 0.5 mm）⇒ **把跖面压到接触**
2. **Step 2**（原 step）：**去掉** prescribe（step 局部 ⇒ 自动失效）、**加**原力载荷
   ⇒ 「**replace the displacement with a force**」

`<Boundary>`（pad 固定）与 `<Contact>` 保持 **global**（两步都生效）。

用法::

    .venv\\Scripts\\python.exe scripts\\opensim_fe\\two_step_pretouch_probe.py --dump
    .venv\\Scripts\\python\\python.exe scripts\\opensim_fe\\two_step_pretouch_probe.py --run
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from pathlib import Path

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

DECK_DIR = ROOT / "temp" / "opensim_fe" / "plantar_fix" / "g0_baseline"
OUT_DIR = ROOT / "temp" / "opensim_fe" / "two_step"
#: 预压位移（mm，沿 −y）。须 > 初始间隙 0.5 mm。
PRETOUCH_MM = 1.0
#: 预压步的时间步数（更细的 ramp ⇒ 更稳）。
PRETOUCH_STEPS = 10
#: `<min_residual>`（手册 §8.3.2.1）。**实测 FEBio 4.13 不认此标签**
#: （顶层/step 级 Control 都报 `unrecognized tag`）⇒ 默认 `None`＝不写。
#: 留作记录：手册给出的「零残差卡住」处方在本版本不可用。
MIN_RESIDUAL = 1e-15
#: **诊断开关**：`True` ⇒ `<Loads>` 留在**全局**（两步都生效）。
#: 用途＝单变量判定「Step 2 残差 ~1e-16」到底是**载荷没进求解**还是别的。
LOADS_GLOBAL = False
#: **诊断开关**：`True` ⇒ 删掉**顶层 `<Control>`**（手册 §6.1：多步分析里 Control
#: 不再定义在文件顶部；pyfebio 的 deck 却有，是唯一剩下的结构偏离）。
DROP_ROOT_CONTROL = False
#: **诊断开关**：`True` ⇒ Step 2 的 `<Control>` 换成**最简形式**
#: （`<time_steps>`＋`<solver/>`）—— 最小复现器 `min_two_step_load_test.py` 用的就是它，
#: 且**载荷确实生效**（Step 2 初始残差 38.9）。
MINIMAL_CTRL = False
#: Step 2 的时间步数。
FORCE_STEPS = 1


def _dump(t: str) -> None:
    print('bytes', len(t))
    for tag in ('<Boundary>', '</Boundary>', '<Loads>', '</Loads>', '<Contact>',
                '</Contact>', '<Steps>', '</Steps>', '<LoadData>', '</LoadData>'):
        idx = [m.start() for m in re.finditer(re.escape(tag), t)]
        print('  %-12s %s' % (tag, idx[:5]))
    i = t.find('<Steps>')
    print('\n--- <Steps> 区段 ---')
    print(t[i:i + 1200])


def _bc(dof: str, val: float) -> str:
    """一条 `prescribed displacement` BC（FEBio 手册 §3.11.6 语法）。"""
    return (
        '        <bc name="pretouch_%s" type="prescribed displacement" '
        'node_set="subtalar_joint">\r\n'
        '          <dof>%s</dof>\r\n'
        '          <value lc="1">%r</value>\r\n'
        '          <relative>0</relative>\r\n'
        '        </bc>\r\n' % (dof, dof, float(val))
    )


def _ctrl_with(ctrl: str, n_steps: int, min_residual: float | None = None) -> str:
    """改写 `<Control>`：``<time_steps>`` 与（可选）``<min_residual>``。

    **`min_residual` 的依据（FEBio 手册 §8.3.2.1「No loads applied」）**：

    > …when the *residual norm* is smaller than a **user-defined limit** FEBio assumes that no
    > loads are present and moves on to the next time step. However, this limit is actually
    > problem dependent and it can happen that for your particular problem the limit is set too
    > small. In that case, the solution is to **increase the limit**.

    本探针实测：预压步成功后 Step 2 的残差 **~1e-16**（实质为零），而 deck 的
    ``min_residual`` 是 **1e-20** ⇒ 卡住不收敛（这不是物理参数，是**零残差判定的数值口径**）。
    """
    out = re.sub(r'<time_steps>\s*\d+\s*</time_steps>',
                 '<time_steps>%d</time_steps>' % int(n_steps), ctrl, count=1)
    if min_residual is not None:
        # ★ 关键：`min_residual` **不是 `<Control>` 的子标签**（FEBio 4.13 报
        # `unrecognized tag`），而是 **`<solver>` 块的参数** —— 依据＝FEBio 运行日志把它
        # 印在 **`SOLVER PARAMETERS`** 标题下；SDK 侧对应 `FENewtonSolver::m_Rmin`
        # （`FECore/FENewtonSolver.h:229`）。手册 §8.3.2.1 说「Control section」是**误导**。
        # deck 的 Control 里正好是空的 `<solver/>` ⇒ 展开它。
        blk = out
        if re.search(r'<solver\s*/>', blk):
            blk = re.sub(r'<solver\s*/>',
                         '<solver>\r\n          <min_residual>%g</min_residual>\r\n'
                         '        </solver>' % float(min_residual), blk, count=1)
        elif re.search(r'<solver>', blk):
            blk = re.sub(r'(<solver>)',
                         r'\1\r\n          <min_residual>%g</min_residual>'
                         % float(min_residual), blk, count=1)
        else:
            raise AssertionError('步内 Control 找不到 <solver>')
        out = blk
    return out


def _root_ctrl(t: str, min_residual: float | None) -> str:
    """在**顶层 `<Control>`** 里插入/改写 `<min_residual>`（手册 §8.3.2.1）。

    ⚠️ **实测 FEBio 4.13 不认此标签** ⇒ `min_residual=None` 时本函数为**恒等**。
    """
    if min_residual is None:
        return t
    i = t.find('<Control>')
    assert i > 0, '找不到顶层 <Control>'
    j = t.find('</Control>', i)
    assert j > i, '找不到顶层 </Control>'
    blk = t[i:j]
    if re.search(r'<min_residual>', blk):
        blk2 = re.sub(r'<min_residual>\s*[^<]*</min_residual>',
                      '<min_residual>%g</min_residual>' % float(min_residual),
                      blk, count=1)
    else:
        blk2 = blk + '    <min_residual>%g</min_residual>\r\n  ' % float(min_residual)
    return t[:i] + blk2 + t[j:]


def _drop_root_control(t: str) -> str:
    """删除**顶层 `<Control>`**（手册 §6.1：多步分析里它不该存在）。"""
    i = t.find('<Control>')
    assert i > 0, '找不到顶层 <Control>'
    j = t.find('</Control>', i)
    assert j > i, '找不到顶层 </Control>'
    j += len('</Control>')
    while j < len(t) and t[j] in '\r\n':
        j += 1
    return t[:i] + t[j:]


def _surgery(t: str) -> str:
    """单步 deck → 两步 deck（手册 §6.1 / §6.1.1 / §6.1.2）。

    原 deck 形态::

        <MeshDomains>…  <Step>  <step id="1" name="Step">
            <Control>…</Control>
            <Boundary>…</Boundary>     <- pad 固定
            <Loads>…</Loads>           <- 压力载荷
            <Contact>…</Contact>
        </step>  </Step>  <Output>…</Output>

    目标形态（§6.1.1：`<Step>` 之前 = **global**，各步都生效；step 内 = 只对该步生效）::

        <Boundary>…</Boundary>          <- pad 固定（global，两步都生效）
        <Contact>…</Contact>            <- 接触（global）
        <Steps>
          <step id="1" name="pretouch">
            <Control>…</Control>
            <Boundary><prescribe bc="y">-P</prescribe></Boundary>   <- 预压到接触
          </step>
          <step id="2" name="force">
            <Control>…</Control>
            <Loads>…</Loads>            <- 「replace the displacement with a force」
          </step>
        </Steps>
        <Output>…</Output>
    """
    m_head = re.search(r'\r?\n\s*<Step>\s*\r?\n', t)
    assert m_head, '找不到 <Step>'
    head = t[:m_head.start()]

    m_tail = re.search(r'</Step>', t)
    assert m_tail, '找不到 </Step>'
    tail = t[m_tail.end():]

    def _block(tag: str, after: int = 0) -> str:
        m = re.search(r'<%s>[\s\S]*?</%s>' % (tag, tag), t[after:])
        assert m, '找不到 <%s>' % tag
        return m.group(0)

    # ⚠️ 必须取 **`<step>` 之内** 的 `<Control>`：deck 顶层另有一个 `<Control>`
    # （含 FEBio 4.13 拒绝的 `<time_stepper type=...>`）⇒ 用最前面的会误抓。
    i_step = t.find('<step ')
    assert i_step > 0, '找不到 <step>'
    ctrl = _block('Control', after=i_step)
    bc = _block('Boundary')          # pad 固定 -> global
    contact = _block('Contact')      # -> global
    loads = _block('Loads')          # -> 只进 step 2

    ind = '      '

    def _ind(blk: str, extra: str = '  ') -> str:
        return (extra + blk.replace('\r\n', '\r\n' + extra)).rstrip()

    step1 = (
        '    <step id="1" name="pretouch">\r\n'
        + _ind(_ctrl_with(_block('Control', after=i_step), PRETOUCH_STEPS))
        + '\r\n      <Boundary>\r\n'
        # ⚠️ 正确标签是 §3.11.6 的 `prescribed displacement`（**不是** `<prescribe>`，
        # FEBio 4.13 报 `unrecognized tag`）；子标签为 `<dof>/<value lc=…>/<relative>`。
        # 语法二次确认自本仓库 `scripts/opensim_fe/g4_fe_stab.py:221`。
        #
        # ⚠️⚠️ **必须三向都控**：只控 `y` 时 x/z 的**刚体模态仍在** ⇒ 依然欠约束
        # （实测：仅 `y` ⇒ 14–44 negJ 全失败）。§8.5.3 的目标是让问题**良约束**，
        # 所以预压步把 `subtalar_joint` 面**完全 prescribed**（x/z = 0、y = −P），
        # Step 2 该面完全释放（step 局部 BC 自动失效）⇒ 回到纯力控。
        + _bc('x', 0.0)
        + _bc('y', -PRETOUCH_MM)
        + _bc('z', 0.0)
        + '      </Boundary>\r\n'
        + '    </step>\r\n'
    )
    _c2 = _block('Control', after=i_step)
    if MINIMAL_CTRL:
        _c2 = ('<Control>\r\n        <time_steps>%d</time_steps>\r\n'
               '        <solver/>\r\n      </Control>' % int(FORCE_STEPS))
    step2 = (
        '    <step id="2" name="force">\r\n'
        + _ind(_ctrl_with(_c2, FORCE_STEPS, min_residual=MIN_RESIDUAL))
        + '\r\n'
        + ('' if LOADS_GLOBAL else _ind(loads) + '\r\n')
        + '    </step>\r\n'
    )
    # 实测：FEBio 4.13 **不认 `<Steps>`**（`unrecognized tag`）；且 `<Step>` 容器内
    # **必须**以 `<step>` 开头（`missing tag "step"`）⇒ 全局段放**顶层**、容器用 `<Step>`。
    # `<value lc="1">` 需要一条 load curve ⇒ 加 `<LoadData>`（顶层，逐字取
    # `scripts/opensim_fe/g4_fe_stab.py` 的 `_RAMP_XML`）。
    ramp = (
        '  <LoadData>\r\n'
        '    <load_controller id="1" type="loadcurve">\r\n'
        '      <interpolate>LINEAR</interpolate>\r\n'
        '      <extend>CONSTANT</extend>\r\n'
        '      <points><point>0,0</point><point>1,1</point></points>\r\n'
        '    </load_controller>\r\n'
        '  </LoadData>\r\n'
    )
    if DROP_ROOT_CONTROL:
        head = _drop_root_control(head)
    return (head
            + '\r\n  ' + bc + '\r\n'
            + '\r\n  ' + contact + '\r\n'
            + (('\r\n  ' + loads + '\r\n') if LOADS_GLOBAL else '')
            + '\r\n' + ramp
            + '\r\n  <Step>\r\n' + step1 + step2 + '  </Step>'
            + tail)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--dump', action='store_true')
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--deck', type=Path, default=DECK_DIR / 'm.feb')
    ap.add_argument('--pretouch-mm', type=float, default=PRETOUCH_MM)
    ap.add_argument('--pretouch-steps', type=int, default=PRETOUCH_STEPS)
    ap.add_argument('--min-residual', type=float, default=MIN_RESIDUAL)
    ap.add_argument('--loads-global', action='store_true')
    ap.add_argument('--drop-root-control', action='store_true')
    ap.add_argument('--minimal-ctrl', action='store_true')
    ap.add_argument('--out-name', type=str, default='m_2step.feb')
    args = ap.parse_args(argv)
    globals()['PRETOUCH_MM'] = float(args.pretouch_mm)
    globals()['PRETOUCH_STEPS'] = int(args.pretouch_steps)
    globals()['MIN_RESIDUAL'] = (None if args.min_residual is None
                                 else float(args.min_residual))
    globals()['LOADS_GLOBAL'] = bool(args.loads_global)
    globals()['DROP_ROOT_CONTROL'] = bool(args.drop_root_control)
    globals()['MINIMAL_CTRL'] = bool(args.minimal_ctrl)

    t = Path(args.deck).read_bytes().decode('iso-8859-1', errors='replace')
    if args.dump:
        _dump(t)
        print('\n=== 手术后 ===')
        print(_surgery(t)[:2400])

    if args.run:
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        out = OUT_DIR / args.out_name
        for s in ('.xplt', '.log', '.hdf5'):
            p = out.with_suffix(s)
            if p.exists():
                p.unlink()
        out.write_bytes(_surgery(t).encode('iso-8859-1'))
        print('[deck]', out)
        # 闸门 3：受控量自证（直接读 deck）
        d = out.read_bytes().decode('iso-8859-1', errors='replace')
        print('  [闸门3] steps=%d  prescribe=%d  Contact=%d  Loads=%d' % (
            len(re.findall(r'<step\b', d)), len(re.findall(r'<prescribe', d)),
            len(re.findall(r'<Contact\b', d, re.I)),
            len(re.findall(r'<Loads\b', d, re.I))))

        from climbing.coupling.febio_run import find_febio
        exe = find_febio()
        t0 = time.time()
        proc = subprocess.run([str(exe), '-i', str(out)], capture_output=True,
                              text=True, cwd=str(OUT_DIR), timeout=1800)
        el = time.time() - t0
        log = out.with_suffix('.log')
        txt = log.read_text(encoding='utf-8', errors='replace') if log.is_file() else ''
        nj = re.findall(r'(\d+)\s*negative\s+jacobians?', txt, re.I)
        print('  [solve] rc=%d  %.1fs  negJ=%s  xplt=%dB' % (
            proc.returncode, el, nj[-1] if nj else None,
            out.with_suffix('.xplt').stat().st_size if out.with_suffix('.xplt').is_file() else 0))
        (OUT_DIR / 'result.json').write_text(json.dumps({
            'pretouch_mm': PRETOUCH_MM, 'rc': proc.returncode,
            'elapsed_s': round(el, 1), 'negJ': nj[-1] if nj else None,
        }, ensure_ascii=False, indent=2), encoding='utf-8')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
