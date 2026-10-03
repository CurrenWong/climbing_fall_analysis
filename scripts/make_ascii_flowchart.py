"""把仿真流程渲染成 ASCII 图，并注入 docs/方案更新_v2.md。

为什么用脚本生成而不是手画
--------------------------
中文字符在等宽字体里占 **2 列**，ASCII 占 1 列。手画的框线一旦有一个字
数错，整张图就错位，而且**改一处要重排全图**。

这里做成生成器：
* ``display_width()`` 按 East Asian Width 计算显示宽度
  （``W``/``F`` 记 2 列，其余记 1 列 —— 框线字符 U+2500 段属
  Ambiguous，按 1 列处理，与"西文等宽 + CJK 回退"的常见渲染一致）
* 画布宽度由内容**自动推导**，改内容不用手动重排
* 只允许 ASCII + 全角字符进入内容（``check_chars``），
  避免混入 »（``·`` ``×`` ``→`` ``①`` ``ε`` 这些 Ambiguous 字符）
  在部分字体下突然变 2 列而错位

注入用标记块，可重复运行（幂等）：
    <!-- ASCII-FLOWCHART:BEGIN -->
    ...
    <!-- ASCII-FLOWCHART:END -->

`tests/test_climbing.py::TestAsciiFlowchart` 会校验文档里的块与生成结果一致，
**防止图与代码漂移**。
"""

from __future__ import annotations

import pathlib
import unicodedata

ROOT = pathlib.Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "方案更新_v2.md"
BEGIN = "<!-- ASCII-FLOWCHART:BEGIN -->"
END = "<!-- ASCII-FLOWCHART:END -->"

#: 两列布局之间的空隙
GAP = 6

# 框线字符
TL, TR, BL, BR = "┌", "┐", "└", "┘"
H, V = "─", "│"
LT, RT = "├", "┤"


# --------------------------------------------------------------------------
# 宽度与字符校验
# --------------------------------------------------------------------------
def display_width(s: str) -> int:
    """按 East Asian Width 计算显示列数（W/F -> 2，其余 -> 1）。"""
    return sum(2 if unicodedata.east_asian_width(c) in ("W", "F") else 1
               for c in s)


def check_chars(s: str, where: str = "") -> None:
    """内容里只允许 ASCII 与全角字符（W/F）。

    混进 Ambiguous 字符（``·`` ``×`` ``→`` ``①②`` 希腊字母 ``ε σ ξ κ``）
    会在部分等宽字体里变成 2 列 —— 图就错位了，且**只在那台机器上错位**。
    """
    bad = [c for c in s
           if not c.isascii() and unicodedata.east_asian_width(c) not in ("W", "F")]
    if bad:
        raise ValueError(
            f"{where}: 含 Ambiguous 宽度字符 {sorted(set(bad))!r} —— "
            f"请改用 ASCII 或中文（如 -> 代替 →、* 代替 ·）"
        )


def pad(s: str, width: int) -> str:
    n = width - display_width(s)
    if n < 0:
        raise ValueError(f"内容超出宽度 {width}（当前 {display_width(s)}）: {s!r}")
    return s + " " * n


# --------------------------------------------------------------------------
# 单个框
# --------------------------------------------------------------------------
def _inner_width(title: str, lines: list[str]) -> int:
    need_title = display_width(title) + 3      # "┌─ title ─…┐"
    need_lines = max((display_width(l) for l in lines), default=0)
    return max(need_title, need_lines) + 2     # 左右各留 1 列内边距


def box(title: str, lines: list[str], inner: int) -> list[str]:
    """渲染一个框，``inner`` 为两竖线之间的列数。"""
    check_chars(title, "title")
    for l in lines:
        check_chars(l, "line")
    k = inner - 3 - display_width(title)
    if k < 0:
        raise ValueError(f"标题过长: {title!r} (inner={inner})")
    out = [f"{TL}{H} {title} " + H * k + TR]
    for l in lines:
        out.append(f"{V}" + pad(l, inner) + V)
    out.append(BL + H * inner + BR)
    # 局部不变量：框内每行等宽，否则右边框会参差
    for l in out:
        if display_width(l) != inner + 2:
            raise ValueError(f"框 {title!r} 内行宽不一致: {l!r}")
    return out


def note(lines: list[str], inner: int) -> list[str]:
    """无标题的框（用于旁注）。"""
    for l in lines:
        check_chars(l, "note")
    out = [LT + H * inner + RT]
    for l in lines:
        out.append(V + pad(l, inner) + V)
    out.append(BL + H * inner + BR)
    for l in out:
        if display_width(l) != inner + 2:
            raise ValueError(f"note 框内行宽不一致: {l!r}")
    return out


def _pad_box(lines: list[str], inner: int, rows: int) -> list[str]:
    """把框补齐到 ``rows`` 行 —— 在底边**之前**插入空行。

    并排的两个框内容行数往往不同（例如「绳索判读」3 行、「抱石判读」4 行）。
    直接 zip 会 IndexError，必须先把矮的补到一样高。
    """
    if len(lines) >= rows:
        return lines
    blank = [V + " " * inner + V] * (rows - len(lines))
    return lines[:-1] + blank + lines[-1:]


def merge_side_by_side(left: list[str], right: list[str], inner: int,
                       gap: int = GAP) -> list[str]:
    rows = max(len(left), len(right))
    left = _pad_box(left, inner, rows)
    right = _pad_box(right, inner, rows)
    return [left[i] + " " * gap + right[i] for i in range(rows)]


# 左右并排时各自的内宽
def _half_inner(pair) -> int:
    a, b = pair
    return max(_inner_width(a[0], a[1]), _inner_width(b[0], b[1]))


# --------------------------------------------------------------------------
# 流程图内容
# --------------------------------------------------------------------------
FLOW: list[tuple] = [
    ("box", "输入参数", [
        "  绳索: 质量 m / 绳长 L0 / 坠落系数 FF / 保护器库",
        "  抱石: 落差 h / 质量 m / 垫参数 / 姿势库 POSTURES",
    ]),

    ("down", None),

    ("box", "人体模型   (B-2 双质点, 主用; src/climbing/pad.py)", [
        "  双质点  m_low (第一接触段) / m_up (躯干)",
        "  m_low_frac 由人体测量学约束, 守卫在 (0,1) 开区间",
        "  屈曲    踝+膝+髋 串联等效弹簧 k_flex + 行程 x_max",
        "  接触    面积 A(塌陷进度) 随下陷时变展开, 不是定值",
        "  姿势    10 种 = 5 个参数, 共用同一套 ODE (加姿势只加数据)",
        "  另有 B-3 六节点链 bodies.py: 结构完成、未标定、当前不用",
    ]),

    # 注意：这里不要再插 ("down", None) —— 紧随其后的 ("split", None)
    # 自己就会发一条竖线 + 分叉线。多一条会渲染出重复的 │ 与 ▼。
    ("split", None),

    ("cols", [
        ("绳索建模   rope.py", [
            "  状态  y = [x, v, L_free]",
            "  几何约束  e = x - h_free",
            "            + (L_free - L0)",
            "  绳本构  T = T_ref * r^n",
            "          r = 应变/应变_ref",
            "  保护器两态: 锁死 / 滑动",
            "             (滑动时限幅 t_hold)",
            "  阻尼  c * de/dt   (仅拉伸耗能)",
            "  运动方程  dv/dt = g - T/m",
        ]),
        ("抱石建模   pad.py", [
            "  状态  y = [z_f, v_f, z_t, v_t]",
            "  双质点  m_low / m_up",
            "          (split_mass 守卫)",
            "  屈曲  xf = L0 - (z_f - z_t)",
            "  垫本构  平台 + 压实",
            "          + 有限压穿支撑",
            "  渐进接触  A(塌陷进度)",
            "            必须时变",
            "  力的配平  见下",
        ]),
    ]),

    # 注意左右：`力的配平` 是**抱石**的概念（被上面 `抱石建模` 的
    # "力的配平 见下" 引用），必须放在**右**列。原 flowchart.png 把它
    # 放在了左列（绳索侧）—— 那是一处内容错位，这里纠正。
    # `冲击窗口定义` 是两模型共用，放左列。
    ("cols", [
        ("冲击窗口定义 (两模型共用)", [
            "  首次接触 t0 -> 首次分离 t1",
            "  峰值只在窗口内取 max",
            "  避免回弹重绷污染",
        ]),
        ("力的配平   (抱石)", [
            "  az_f = g - f_pad/m_low",
            "          + f_flex/m_low",
            "  az_t = g - f_flex/m_up",
            "  两式相加 = f_pad/m - g",
            "  (质心方程, 牛顿第三定律)",
        ]),
    ]),

    ("join", "组装 ODE"),

    ("box", "数值求解   solve_ivp", [
        "  绳索: RK45   max_step 2e-4   rtol 1e-8",
        "  抱石: RK45   max_step 1e-4   rtol 1e-8   atol 1e-11",
        "  (B-3 链模型改用隐式 Radau -- 刚性接触)",
    ]),

    ("down", None),

    ("box", "后处理  (窗口内统一处理)", [
        "  峰值    max(力) / max(加速度) / max(接触压力 p = F/A)",
        "          全部在窗口内取, 不是全程",
        "  能量账  用功率累积的峰值",
        "          回弹段会把能量减回来; trapezoid(f, x)",
        "          对弹性体净积分为 0 -> 算出负数",
        "          E_pad 与 E_flex 取同一时刻快照",
        "          相加才不会超过输入能量",
        "  单位    accel 统一除 g; 字段名以 _g 结尾",
        "          就必须真的是 g",
    ]),

    ("down", None),

    ("cols", [
        ("绳索指标", [
            "  峰值张力 / 加速度 / 伸长率",
            "  放绳量 / 绳总滑动 / 滑动态时长",
            "  保护者侧张力 = exp(mu*theta)",
            "  即欧拉绞盘放大",
        ]),
        ("抱石指标", [
            "  峰值垫反力 / 接触压力 F/A",
            "  各段加速度 (第一接触部位)",
            "  屈曲行程 / 能量分配 / 残余",
        ]),
    ]),

    ("down", None),

    ("cols", [
        ("绳索判读", [
            "  UIAA 101: 峰值 < 12 kN",
            "  峰值对 FF(>1) 严格单调",
            "  保护器峰值 = min(绳, t_hold)",
        ]),
        ("抱石判读 (分部位)", [
            "  接触部位  压力 150/400/800 kPa",
            "  头/颈     HIC 240/700/1000",
            "  脊柱/下肢  骨载荷 6/10/16 kN",
            "  (v2 P0 已换成 7 个骨部位阈值)",
        ]),
    ]),

    ("labels", ("标定", "验证")),

    ("cols", [
        ("标定   calibrate_rope.py", [
            "  UIAA 条件网格搜索",
            "  目标 = 加权 RSS",
            "         (峰值 .7 / 伸长 .3)",
            "  相对误差相加会正负抵消成假优解",
            "  标定完必须写回默认值",
            "  并检查下游功能是否活过来",
        ]),
        ("验证   tests/  97 用例", [
            "  正交不变量:",
            "    牛顿第三定律 / 能量守恒",
            "    单调性 / 已知解(HIC 闭式)",
            "    / 单位 / 反函数自洽",
            "  约一半用例是锁死已修的静默 bug",
        ]),
    ]),

    ("feedback", "闭环"),

    ("box", "v2 附加模块  (本次方案新增的部分)", [
        "  P0  bone.py      ->  抱石判读: 7 个骨部位阈值",
        "                      (取代笼统的 6/10/16 kN)",
        "  P1  scenarios.py ->  输入参数: 场景权重与采样",
        "                      (取代手写姿势权重)",
    ]),

    ("down", None),

    ("box", "实际跑过的参数矩阵  (明细见 docs/方案更新_v2.md 1.3 节)", [
        "  用途            高度范围       步长    姿势         质量       落点",
        "  B-2 全套分析    0.3 - 5.0 m    9 点    全部 10 种   80 kg      软垫+硬地",
        "  时程/压力图     0.5 - 5.0 m    0.5 m   逐个         80 kg      软垫+硬地",
        "  对标 [Y25] 阈值 1   - 40  m    1 m     绷直腿       77 kg      刚性地面",
        "  P0/P1 蒙特卡洛  0.5 - 4.5 m    随机    B25 权重     40-130 kg  垫94%/硬6%",
        "  绳索 UIAA 101   FF 1.1 - 2.0   网格    -            80 kg      -",
        "  没跑过: 旋转 / 落地角度 / 踝关节 / 上肢 / 滑动 / 垫型号对比",
    ]),
]


# --------------------------------------------------------------------------
# 渲染
# --------------------------------------------------------------------------
def render() -> str:
    single_w = 0
    half_w = 0
    for op in FLOW:
        kind = op[0]
        if kind == "box":
            single_w = max(single_w, _inner_width(op[1], op[2]))
        elif kind == "cols":
            half_w = max(half_w, _half_inner(op[1]))

    # ---- 画布宽度推导 ----
    # 目标恒等式（右边框能对齐**全靠**它们，下方有 assert 守着）：
    #   * 单列框宽度    = single_inner + 2        = box_total
    #   * 右列框右边缘  = right_col + half_total = box_total
    # 第二条要求 (box_total - GAP) 能被 2 整除，故画布按需补到偶数差。
    #
    # 早期版本写的是 half_inner = max(half_w, (box_total-GAP-4)//2)：
    # 当 box_total 被**单列内容**撑大时，half_w 反而成为下限，
    # 于是 half_total 变大 -> right_col + half_total > box_total，
    # **两条恒等式同时失效**而布局照样渲染出来（右边框参差，
    # 看上去仍然"像图"，肉眼极难发现）。现在改成由 GAP 反解。
    half_need = half_w + 2                     # 每侧框至少需要的总宽
    box_total = max(single_w + 2, 2 * half_need + GAP)
    if (box_total - GAP) % 2:                 # 让两列精确均分
        box_total += 1
    half_total = (box_total - GAP) // 2
    half_inner = half_total - 2
    single_inner = box_total - 2
    right_col = half_total + GAP

    spine = box_total // 2          # 单列时的竖线/箭头列
    left_ctr = half_total // 2      # 左侧框的中心列
    right_ctr = right_col + half_total // 2

    def at(col: int, text: str) -> str:
        return " " * col + text

    def two(col_a: int, text_a: str, col_b: int, text_b: str) -> str:
        s = " " * col_a + text_a
        gap = col_b - display_width(s)
        if gap < 0:
            raise ValueError(f"两段文字重叠: {text_a!r} @{col_a} 与 {text_b!r} @{col_b}")
        return s + " " * gap + text_b

    def tee(open_ch: str, close_ch: str, stub: str) -> str:
        """``┌───┴───┐`` 形状的分叉/汇合线。"""
        return (at(left_ctr, open_ch
                   + H * (spine - left_ctr - 1) + stub
                   + H * (right_ctr - spine - 1) + close_ch))

    out: list[str] = []

    for op in FLOW:
        kind = op[0]

        if kind == "box":
            out += box(op[1], op[2], single_inner)

        elif kind == "cols":
            (t1, l1), (t2, l2) = op[1]
            out += merge_side_by_side(box(t1, l1, half_inner),
                                      box(t2, l2, half_inner), half_inner)

        elif kind == "cols_note":
            out += merge_side_by_side(note(op[1][0][0], half_inner),
                                      note(op[1][0][1], half_inner), half_inner)

        elif kind == "split":
            out.append(at(spine, V))
            out.append(tee("┌", "┐", "┴"))
            out.append(two(left_ctr, "▼", right_ctr, "▼"))

        elif kind == "join":
            out.append(tee("└", "┘", "┬"))
            out.append(at(spine - 1, "▼  " + (op[1] or "")))

        elif kind == "down":
            out.append(at(spine, V))
            out.append(at(spine, "▼"))

        elif kind == "labels":
            a, b = op[1]
            out.append(two(left_ctr, V, right_ctr, V))
            out.append(two(left_ctr, "▼ " + a, right_ctr, "▼ " + b))

        elif kind == "feedback":
            out.append(tee("└", "┘", "┬"))
            out.append(at(spine - 1, "▼  " + op[1]))

        else:
            raise ValueError(f"未知操作: {kind}")

    # 布局恒等式 —— 右边框能对齐**靠的是这两条**，在此显式断言。
    #   * 单列框宽度 = single_inner + 2 = box_total
    #   * 右列框右边缘 = right_col + half_total = 2*half_total + GAP = box_total
    assert single_inner + 2 == box_total, (single_inner, box_total)
    assert right_col + half_total == box_total, (right_col, half_total, box_total)
    assert half_inner >= half_w, (half_inner, half_w)
    assert (box_total - GAP) % 2 == 0, box_total

    self_check(out, box_total)
    body = "\n".join(out)
    return ("```text\n"
            "攀岩坠落仿真 · 完整流程\n"
            + body
            + "\n```")


def self_check(lines: list[str], canvas: int) -> None:
    """正交不变量：任何一行都不得超出画布。

    这条检查不是装饰 —— 中文字宽数错一个字，图就错位，
    而**错位的图看上去仍然"像图"**，肉眼很难发现。属于本项目
    `阶段总结.md` §3.1 那类"不会自己暴露"的错误。

    （"右边框对齐"不在这里查：居中的分叉线 ``┌───┴───┐`` 从 left_ctr 起、
    到 right_ctr 止，本来就短于画布。`box()` 内部已断言每行等宽，
    布局恒等式也已在 `render()` 里断言，所以对齐性质是构造保证的。）
    """
    for i, l in enumerate(lines, 1):
        w = display_width(l)
        if w > canvas:
            raise ValueError(f"第 {i} 行超出画布 {canvas}（宽 {w}）: {l!r}")


# --------------------------------------------------------------------------
# 注入文档
# --------------------------------------------------------------------------
def block() -> str:
    return f"{BEGIN}\n{render()}\n{END}"


def inject(doc: pathlib.Path = DOC, dry_run: bool = False) -> str:
    text = doc.read_text(encoding="utf-8")
    new = block()
    if BEGIN in text and END in text:
        head, rest = text.split(BEGIN, 1)
        _, tail = rest.split(END, 1)
        updated = head + new + tail
    else:
        raise RuntimeError(
            f"{doc.name} 里没有 {BEGIN} / {END} 标记块 —— 请先手工放一对标记")
    if not dry_run and updated != text:
        doc.write_text(updated, encoding="utf-8")
    return updated


if __name__ == "__main__":
    import sys
    if "--widths" in sys.argv:
        sw = hw = 0
        for op in FLOW:
            if op[0] == "box":
                w = _inner_width(op[1], op[2])
                sw = max(sw, w)
                print(f"  single {w:3d}  {op[1]}")
                for l in op[2]:
                    if display_width(l) + 2 >= w:
                        print(f"         {display_width(l):3d}  {l!r}")
            elif op[0] == "cols":
                w = _half_inner(op[1])
                hw = max(hw, w)
                print(f"  half   {w:3d}  {op[1][0][0]} | {op[1][1][0]}")
                for t, ls in op[1]:
                    for l in ls:
                        if display_width(l) + 2 >= w:
                            print(f"         {display_width(l):3d}  {l!r}")
        print(f"single_w={sw} half_w={hw}")
        raise SystemExit(0)

    dry = "--dry" in sys.argv
    txt = render()
    print(txt)
    print(f"\n--- 画布宽度: {max(display_width(l) for l in txt.splitlines())} 列 ---")
    if not dry:
        inject()
        print(f"已注入 {DOC.relative_to(ROOT)}")
