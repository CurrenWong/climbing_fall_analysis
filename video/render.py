# -*- coding: utf-8 -*-
"""
《抱石怎么摔才不伤？》短视频渲染器
路线 A：纯数据动画 + 简笔示意 + TTS 旁白 + 硬字幕（无 BGM、无 AI 镜头）

用法:
  python render.py preview   -> 只出 6 张每镜中点帧 PNG（校对视觉用）
  python render.py full      -> 渲染全片并合成音轨/字幕，输出 mp4
"""
import io
import os
import subprocess
import sys
import wave

import matplotlib
matplotlib.use("Agg")
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle, Wedge, FancyArrowPatch
import matplotlib.patheffects as pe
import numpy as np

ROOT = r"D:\Project\climbing_fall_analysis\video"
VO_DIR = os.path.join(ROOT, "vo")
TMP = os.path.join(ROOT, "temp")
OUT = os.path.join(ROOT, "out")
FONT_PATH = r"C:\Windows\Fonts\msyh.ttc"
FONT_BOLD = r"C:\Windows\Fonts\msyhbd.ttc"
FFMPEG = r"D:\Program\ffmpeg-8.0-full_build\bin\ffmpeg.exe"

for p in (FONT_PATH, FONT_BOLD):
    if os.path.exists(p):
        font_manager.fontManager.addfont(p)
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei"]
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["axes.unicode_minus"] = False

W, H, FPS = 1080, 1920, 30
GAP, TAIL = 0.40, 1.50
FS = 3.3  # 字号放大系数：figsize 以英寸计，1080x1920@100dpi 下 pt 需整体放大

# 变体开关（默认行为不变）：
#   S1_VARIANT=stick    第一镜用简笔坠落剪影（默认）
#   S1_VARIANT=kenburns 第一镜改用 wanx 图做 Ken Burns 慢推
#   RENDER_TAG=_kb      产物文件名后缀，避免覆盖已有成片
S1_VARIANT = os.environ.get("S1_VARIANT", "stick")
RENDER_TAG = os.environ.get("RENDER_TAG", "")
KB_IMAGE = r"D:\Project\climbing_fall_analysis\paper\wechat_assets\wanx_raw_1.png"

BG = "#0B0E14"
PANEL = "#161B26"
PANEL2 = "#1E2532"
FG = "#F0F3F8"
MUTED = "#8B96A8"
BLUE = "#4C9AFF"
ORANGE = "#FF7B45"
GREEN = "#3FB950"
RED = "#F85149"
YELLOW = "#E3B341"
PURPLE = "#A371F7"

VO_TEXT = {
    1: "抱石的人，迟早要摔一次。可你真的知道自己是怎么摔的吗？",
    2: "一项法国研究追踪了245名抱石者的301起损伤：67%伤在下肢，踝关节占四成，"
       "伤后不能攀爬的中位时间是90天。",
    3: "更关键的是，坠落的方式分两条完全不同的路。第一条：身体保持垂直、坠落中旋转、"
       "最后脚先着地——占大多数，伤在脚踝。脚和踝陷进垫子，踝被迫过度外旋，韧带首当其冲。"
       "第二条：起手就是后仰，坠落中大幅转体，落地不稳时本能用手撑垫子——伤在肘和肩。",
    4: "最常见的前三种场景里，两种都是垂直落体、脚先着地，加起来超过三成，结果都是下肢伤。",
    5: "什么让人伤得更重？28到31岁的重伤风险，是32岁以上人群的3.35倍。"
       "从墙底摔，比从墙顶摔低82%。",
    6: "所以：双手交叉护住身前，落地顺势后仰翻滚。完整解读，看公众号。",
}


# ---------------------------------------------------------------- helpers
def smooth(x):
    x = max(0.0, min(1.0, x))
    return x * x * (3 - 2 * x)


def seg(t, t0, t1):
    """t 在 [t0,t1] 区间内的平滑进度，区间外返回 0 或 1。"""
    if t1 <= t0:
        return 1.0
    return smooth((t - t0) / (t1 - t0))


def new_fig():
    fig = plt.figure(figsize=(W / 100, H / 100), dpi=100)
    fig.patch.set_facecolor(BG)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 108)
    ax.set_ylim(0, 192)
    ax.axis("off")
    ax.set_facecolor(BG)
    return fig, ax


def txt(ax, x, y, s, size=6, color=FG, weight="bold", ha="center", va="center",
        alpha=1.0, stroke=0, family="Microsoft YaHei", zorder=10, ls="normal"):
    t = ax.text(x, y, s, fontsize=size * FS, color=color, fontweight=weight, ha=ha, va=va,
                alpha=alpha, zorder=zorder, fontfamily=family, linespacing=1.5)
    if ls != "normal":
        t.set_linestyle(ls)
    if stroke:
        t.set_path_effects([pe.withStroke(linewidth=stroke, foreground="#000000")])
    return t


def panel(ax, x, y, w, h, fc=PANEL, ec=None, lw=1.2, r=2.2, alpha=1.0, zorder=2):
    p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                       facecolor=fc, edgecolor=ec or fc, linewidth=lw,
                       alpha=alpha, zorder=zorder)
    ax.add_patch(p)
    return p


def arrow(ax, x0, y0, x1, y1, color=BLUE, lw=4, alpha=1.0, zorder=6, style="-|>"):
    a = FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style,
                        mutation_scale=26, lw=lw, color=color, alpha=alpha, zorder=zorder,
                        shrinkA=0, shrinkB=0)
    ax.add_patch(a)
    return a


def progress_bar(ax, t, dur, color=BLUE):
    """底部细进度条。"""
    y = 2.0
    ax.add_patch(Rectangle((0, y), 108, 0.7, facecolor="#1A2030", edgecolor="none", zorder=20))
    ax.add_patch(Rectangle((0, y), 108 * min(1.0, t / dur), 0.7, facecolor=color,
                           edgecolor="none", zorder=21))


# ---------------------------------------------------------------- subtitles
def wrap(text, width=14):
    """按标点优先换行；仍然超长的行再硬切。"""
    lines, cur = [], ""
    for ch in text:
        cur += ch
        if len(cur) >= width and ch in "，、。：；？！":
            lines.append(cur)
            cur = ""
    if cur:
        lines.append(cur)
    out = []
    for l in lines:
        while len(l) > width:
            out.append(l[:width])
            l = l[width:]
        if l:
            out.append(l)
    return out or [""]


def split_cues(text, max_lines=2, width=14):
    """把旁白切成字幕块：先按句号类标点分句，再按 width 换行，每块最多 max_lines 行。"""
    parts, buf = [], ""
    for ch in text:
        buf += ch
        if ch in "。？！：":
            parts.append(buf)
            buf = ""
    if buf:
        parts.append(buf)
    lines = []
    for p in parts:
        lines += wrap(p, width)
    return ["\n".join(lines[i:i + max_lines]) for i in range(0, len(lines), max_lines)]


def cue_times(cues, dur, lead=0.10):
    """按字数比例把音轨时长分配给各条字幕。"""
    total = sum(len(c) for c in cues) or 1
    out, t = [], lead
    span = dur - lead
    for c in cues:
        d = span * len(c) / total
        out.append((t, t + d, c))
        t += d
    return out


def draw_subtitle(ax, cues, t, dur):
    body = None
    for c0, c1, text in cues:
        if c0 <= t < c1:
            body = text
            break
    txt(ax, 54, 22, body if body else "", size=5.6, color="#FFFFFF", stroke=5)
    progress_bar(ax, t, dur)


# ---------------------------------------------------------------- shots
def shot1(ax, t, dur):
    """钩子：坠落剪影 + 大字钩子"""
    # 背景竖向渐变（深蓝灰，不用 magma 的紫红）
    grad = np.linspace(0, 1, 256).reshape(-1, 1) / 255.0
    rgb = np.zeros((256, 1, 3))
    rgb[:, 0, 0] = 0.42 + 0.22 * grad[:, 0]
    rgb[:, 0, 1] = 0.50 + 0.28 * grad[:, 0]
    rgb[:, 0, 2] = 0.72 + 0.28 * grad[:, 0]
    ax.imshow(rgb, extent=[0, 108, 0, 192], aspect="auto", zorder=0, alpha=0.22)

    # 墙面
    ax.add_patch(Rectangle((0, 0), 16, 192, facecolor="#151A24", edgecolor="none", zorder=1))
    for i in range(16):
        ax.add_patch(Rectangle((2 + (i % 5) * 2.8, 10 + i * 12), 1.6, 1.6,
                               facecolor=ORANGE, edgecolor="none", alpha=0.45, zorder=2))

    # 垫子
    panel(ax, 14, 30, 92, 8, fc="#1B3A5C", ec="#2D5F8F", r=1.6, zorder=3)
    txt(ax, 60, 34, "保护垫", size=3.2, color="#5C8CB8", weight="normal", zorder=4)

    # 坠落剪影（头 + 躯干 + 四肢）
    fall = seg(t, 0.15, dur - 0.2)
    cx = 30 + 44 * fall
    cy = 150 - 108 * (fall ** 1.7)
    ax.add_patch(Circle((cx, cy + 4), 2.0, facecolor=FG, edgecolor="none", zorder=8))
    ax.plot([cx, cx + 2.4], [cy + 2, cy - 4], color=FG, lw=5, solid_capstyle="round", zorder=8)
    ax.plot([cx + 2.4, cx - 1.0], [cy - 4, cy - 9], color=FG, lw=4, solid_capstyle="round", zorder=8)
    ax.plot([cx + 2.4, cx + 6.0], [cy - 4, cy - 8.5], color=FG, lw=4, solid_capstyle="round", zorder=8)
    ax.plot([cx + 1.0, cx + 5.5], [cy + 1, cy + 3.5], color=FG, lw=3.4, solid_capstyle="round", zorder=8)
    ax.plot([cx + 1.0, cx - 2.5], [cy + 1, cy + 3.0], color=FG, lw=3.4, solid_capstyle="round", zorder=8)

    # 速度线：拖在人物上方
    for i in range(6):
        yy = cy + 12 + i * 5.0
        ln = 9 + 7 * ((i * 7 + t * 40) % 5) / 5
        ax.plot([cx - ln, cx], [yy, yy], color="#FFFFFF", lw=1.8,
                alpha=0.10 + 0.26 * ((i * 3 + int(t * 12)) % 4) / 4, zorder=7)

    # 钩子大字
    a1 = seg(t, 0.3, 1.1)
    a2 = seg(t, 1.5, 2.3)
    txt(ax, 54, 172, "抱石的人，迟早要摔一次", size=7.4, alpha=a1)
    txt(ax, 54, 158, "可你真的知道自己是怎么摔的吗？", size=7.4,
        color=YELLOW, alpha=a2)

    a3 = seg(t, 2.6, 3.4)
    panel(ax, 18, 134, 72, 11, fc=PANEL, ec=PANEL2, alpha=a3, zorder=5)
    txt(ax, 54, 139.5, "245 名抱石者 · 301 起损伤", size=5.0, color=MUTED, alpha=a3)


# ---- 变体 B：S1 改用实拍风格图 + Ken Burns 慢推 ----
# 源图 1600x1280（4:3 横构图）塞进 1080x1920 只能裁切，所以先把源图裁到
# z=1 时那个 9:16 窗口（宽 720px、高正好等于源图高），锚点锁在坠落者身上。
KB_ANCHOR_PX = (1112.0, 640.0)  # 源图坐标：坠落者重心
KB_CROP_W = 720                 # z=1 窗口宽度（源图像素）= H/W 比例下的满幅宽度

_KB_MASTER = None
_KB_ANCHOR = None


def _kb_master():
    """裁剪后的母图 + 锚点在其内的坐标。裁掉的部分永远不会出现在画面里。"""
    global _KB_MASTER, _KB_ANCHOR
    if _KB_MASTER is None:
        im = Image.open(KB_IMAGE).convert("RGB")
        iw, ih = im.size
        cx, cy = KB_ANCHOR_PX
        w = min(KB_CROP_W, iw)
        h = min(int(round(KB_CROP_W * H / W)), ih)
        x0 = int(round(min(max(cx - w / 2.0, 0), iw - w)))
        y0 = int(round(min(max(cy - h / 2.0, 0), ih - h)))
        _KB_MASTER = np.asarray(im.crop((x0, y0, x0 + w, y0 + h)))
        _KB_ANCHOR = (cx - x0, cy - y0)
    return _KB_MASTER, _KB_ANCHOR


def vert_scrim(ax, y_lo, y_hi, a_lo, a_hi, zorder=1):
    """竖向黑色渐隐遮罩：底部/顶部压暗，保证叠字可读。"""
    n = 256
    a = np.linspace(a_lo, a_hi, n)
    rgba = np.zeros((n, 1, 4), dtype=float)
    rgba[:, 0, 3] = a
    ax.imshow(rgba, extent=[0, 108, y_lo, y_hi], aspect="auto",
              zorder=zorder, interpolation="bilinear", origin="lower")


def shot1_kenburns(ax, t, dur):
    """钩子（变体 B）：wanx 图 Ken Burns 慢推 + 钩子大字。

    注意 extent 宽度必须 *乘以* z：轴单位越大，画面被撑得越大、露出的源图
    越少，那才是"推近"；除以 z 是往后退，还会让画面缩到填不满屏。
    """
    im, (au, av) = _kb_master()
    ih, iw = im.shape[:2]

    z = 1.0 + 0.12 * seg(t, 0.0, dur)
    mw, mh = 108.0 * z, 192.0 * z
    x0 = 54.0 - (au / iw) * mw
    y0 = 96.0 - (1.0 - av / ih) * mh
    ax.imshow(im, extent=[x0, x0 + mw, y0, y0 + mh], aspect="auto",
              zorder=0, interpolation="lanczos")

    # 整体压暗一点，再接上/下渐隐
    ax.add_patch(Rectangle((0, 0), 108, 192, facecolor="#000000",
                           edgecolor="none", alpha=0.22, zorder=1))
    vert_scrim(ax, 0, 42, 0.78, 0.0, zorder=1)
    vert_scrim(ax, 130, 192, 0.0, 0.82, zorder=1)

    a1 = seg(t, 0.3, 1.1)
    a2 = seg(t, 1.5, 2.3)
    txt(ax, 54, 172, "抱石的人，迟早要摔一次", size=7.4, alpha=a1, stroke=6)
    txt(ax, 54, 158, "可你真的知道自己是怎么摔的吗？", size=7.4,
        color=YELLOW, alpha=a2, stroke=6)

    a3 = seg(t, 2.6, 3.4)
    panel(ax, 18, 134, 72, 11, fc=PANEL, ec=PANEL2, alpha=a3 * 0.90, zorder=5)
    txt(ax, 54, 139.5, "245 名抱石者 · 301 起损伤", size=5.0, color=MUTED,
        alpha=a3, zorder=6)


def shot2(ax, t, dur):
    """数据卡：伤情分布堆叠条 + 三个数字块 + 严重度"""
    txt(ax, 54, 180, "伤在哪里", size=7.6)
    txt(ax, 54, 170, "301 起损伤的部位分布", size=4.0, color=MUTED, weight="normal")

    # 部位分布（从左到右：大 → 小），游标只前进不回退
    segs = [("下肢", 67, BLUE), ("上肢", 25, ORANGE), ("躯干", 4, MUTED), ("头颈", 3, PURPLE)]
    bx, bw, by, bh = 10, 88, 150, 7
    cursor = bx
    for i, (name, val, col) in enumerate(segs):
        frac = val / 100
        local = seg(t, 0.3 + i * 0.22, 1.1 + i * 0.22)
        ax.add_patch(Rectangle((cursor, by), bw * frac * local, bh,
                               facecolor=col, edgecolor="none", zorder=4))
        if frac >= 0.12 and local > 0.7:
            txt(ax, cursor + bw * frac / 2, by + bh / 2, f"{val}%", size=4.2,
                color="#0B0E14", alpha=seg(t, 0.6 + i * 0.22, 1.3 + i * 0.22), zorder=6)
        cursor += bw * frac
    ax.add_patch(Rectangle((bx, by), bw, bh, facecolor="none", edgecolor="#2A3242",
                           linewidth=1.0, zorder=7))

    # 图例（含两个窄分段，避免标签叠字）
    lx = bx
    la = seg(t, 0.9, 1.7)
    for name, val, col in segs:
        ax.add_patch(Rectangle((lx, 142.5), 2.2, 2.2, facecolor=col,
                               edgecolor="none", alpha=la, zorder=5))
        txt(ax, lx + 3.2, 143.6, f"{name} {val}%", size=3.2, color=MUTED,
            weight="normal", alpha=la, ha="left", zorder=5)
        lx += 22

    # 三个数字块
    tiles = [
        ("踝关节", "40%", "全部损伤中占比", ORANGE),
        ("扭伤", "36%", "最常见类型", BLUE),
        ("停训中位", "90 天", "伤后不能攀爬", RED),
    ]
    tw, th, ty = 27, 34, 104
    for i, (lab, val, note, col) in enumerate(tiles):
        tx = 10 + i * (tw + 3.5)
        loc = seg(t, 2.0 + i * 0.3, 2.8 + i * 0.3)
        panel(ax, tx, ty, tw, th, fc=PANEL, ec=PANEL2, alpha=0.3 + 0.7 * loc, zorder=3)
        txt(ax, tx + tw / 2, ty + th - 8, lab, size=3.8, color=MUTED,
            weight="normal", alpha=loc)
        if lab == "停训中位":
            shown = f"{int(round(90 * loc))} 天"
        else:
            shown = f"{int(round(float(val.rstrip('%')) * loc))}%"
        txt(ax, tx + tw / 2, ty + th / 2 - 1, shown, size=6.6, color=col, alpha=loc)
        txt(ax, tx + tw / 2, ty + 6, note, size=3.2, color=MUTED, weight="normal", alpha=loc)

    a = seg(t, 3.2, 3.9)
    txt(ax, 54, 88, "伤后不能爬的时间中位数，接近三个月", size=4.6, color=YELLOW, alpha=a)

    # 严重度分布
    sev = [("UIAA 1 级", 14, MUTED), ("2 级", 64, YELLOW), ("3 级及以上", 20, RED)]
    b = seg(t, 4.2, 5.0)
    panel(ax, 10, 48, 88, 30, fc=PANEL, ec=PANEL2, alpha=0.3 + 0.7 * b, zorder=3)
    txt(ax, 18, 70, "严重度", size=3.8, color=MUTED, weight="normal", alpha=b, ha="left")
    cur = 18
    for name, val, col in sev:
        w = 72 * val / 98
        ax.add_patch(Rectangle((cur, 58), w * b, 4, facecolor=col,
                               edgecolor="none", zorder=5))
        if w > 10:
            txt(ax, cur + w / 2, 60, f"{val}%", size=3.4, color="#0B0E14",
                alpha=b, zorder=6)
        cur += w
    txt(ax, 18, 52.5, "3 级及以上（重伤）占 20%", size=3.4, color=RED,
        alpha=b, ha="left")


def _path_col(ax, x, title, color, rows, foot, t, dur):
    panel(ax, x, 46, 46, 96, fc=PANEL, ec=color, lw=1.6, alpha=0.95, zorder=3)
    txt(ax, x + 23, 133, title, size=5.0, color=color, zorder=10)
    y = 116
    for i, (lab, val) in enumerate(rows):
        loc = seg(t, 0.6 + i * 1.5, 1.6 + i * 1.5)
        ax.add_patch(Circle((x + 8, y), 3.6, facecolor=color, edgecolor="none",
                            alpha=0.22 + 0.78 * loc, zorder=6))
        txt(ax, x + 8, y, str(i + 1), size=3.4, color="#0B0E14", alpha=loc, zorder=7)
        txt(ax, x + 14, y, lab, size=4.0, ha="left", alpha=loc, zorder=7)
        if val:
            txt(ax, x + 41, y, val, size=5.0, color=color, alpha=loc, ha="right", zorder=7)
        if i < len(rows) - 1:
            ax.plot([x + 8, x + 8], [y - 4.0, y - 20.0], color=color, lw=2.4,
                    alpha=0.7 * loc, zorder=5)
            arrow(ax, x + 8, y - 17.0, x + 8, y - 20.2, color=color, lw=2.4,
                  alpha=0.8 * loc, zorder=6)
        y -= 26
    f = seg(t, dur * 0.62, dur * 0.74)
    ax.plot([x + 4, x + 42], [66, 66], color=color, lw=1.0, alpha=0.4 * f, zorder=5)
    txt(ax, x + 23, 57, foot[0], size=4.2, color=color, alpha=f, zorder=7)
    txt(ax, x + 23, 51, foot[1], size=3.2, color=MUTED, weight="normal", alpha=f, zorder=7)


def shot3(ax, t, dur):
    """核心镜：左右分屏，两条坠落路径"""
    a = seg(t, 0.0, 0.8)
    txt(ax, 54, 178, "两种摔法，两种伤", size=7.6, alpha=a)
    txt(ax, 54, 168, "下降肢 vs 上肢，走的是两条完全不同的路径", size=4.0,
        color=MUTED, weight="normal", alpha=a)

    _path_col(ax, 5, "路径一 · 伤在脚踝", BLUE,
              [("身体保持垂直", "73%"), ("坠落中旋转", "61%"), ("最后脚先着地", "87%")],
              ("踝扭伤占 28%", "脚踝陷进垫子，被迫过度外旋"), t, dur)
    _path_col(ax, 57, "路径二 · 伤在肘和肩", ORANGE,
              [("起手就后仰", ""), ("坠落中大幅转体", ""), ("本能用手撑垫子", "")],
              ("上肢伤占 11%", "其中 5 例肘脱位"), t, dur)

    c = seg(t, dur * 0.80, dur * 0.90)
    panel(ax, 5, 30, 98, 12, fc=PANEL2, ec=PANEL, lw=1.0, alpha=c, zorder=3)
    txt(ax, 54, 36, "同样的一个坠落，落在哪儿决定伤哪儿", size=4.4, alpha=c)


def shot4(ax, t, dur):
    """最高频的两种场景"""
    txt(ax, 54, 180, "最高频的两种摔法", size=7.4)
    txt(ax, 54, 170, "加起来占全部场景的 31%", size=4.0, color=MUTED, weight="normal")

    cards = [
        ("站立着地", "42 例", "17.1%", "90% 落在下肢", BLUE),
        ("脚先着地 + 后仰", "35 例", "14.3%", "77% 落在下肢", BLUE),
    ]
    for i, (name, cnt, pct, note, col) in enumerate(cards):
        y = 112 - i * 50
        loc = seg(t, 0.4 + i * 0.5, 1.4 + i * 0.5)
        panel(ax, 12, y, 84, 44, fc=PANEL, ec=PANEL2, alpha=0.25 + 0.75 * loc, zorder=3)
        txt(ax, 22, y + 33, f"{i+1}", size=5.6, color=col, ha="left", alpha=loc)
        txt(ax, 32, y + 33, name, size=4.8, ha="left", alpha=loc)
        txt(ax, 32, y + 24, f"{cnt} · {note}", size=3.6, color=MUTED,
            weight="normal", ha="left", alpha=loc)
        f = float(pct.rstrip("%")) / 20.0
        w = 62 * f * loc
        ax.add_patch(Rectangle((32, y + 12), w, 5, facecolor=col, edgecolor="none", zorder=4))
        txt(ax, 93, y + 14.5, pct, size=5.6, color=col, ha="right", alpha=loc)

    a = seg(t, 1.8, 2.5)
    panel(ax, 12, 34, 84, 22, fc="#12263F", ec=BLUE, lw=1.4, alpha=a, zorder=3)
    txt(ax, 54, 49, "两种都是垂直落体 · 脚先着地", size=4.4, color=BLUE, alpha=a)
    txt(ax, 54, 40, "结果都是下肢伤 —— 踝是第一重灾区", size=3.6,
        color=FG, alpha=a)


def shot5(ax, t, dur):
    """森林图：谁伤得更重"""
    txt(ax, 54, 178, "谁伤得更重", size=7.4)
    txt(ax, 54, 167, "比值（OR）· 虚线 = 无差异", size=4.6, color=MUTED, weight="normal")

    rows = [
        ("28–31 岁", "3.35", "vs 32 岁以上", RED, 1.0),
        ("中等水平", "0.27", "vs 困难", GREEN, 0.0),
        ("从墙底摔", "0.18", "vs 从墙顶摔", GREEN, 0.0),
    ]
    # log10 刻度 0.1 -> 6
    lo, hi = np.log10(0.1), np.log10(6.0)
    x0, x1, yy0 = 30, 88, 146

    def px(v):
        return x0 + (np.log10(v) - lo) / (hi - lo) * (x1 - x0)

    for gv in (0.1, 0.3, 1.0, 3.0):
        ax.plot([px(gv), px(gv)], [yy0 - 4, yy0 - 46], color="#232B3A", lw=1, zorder=1)
        txt(ax, px(gv), yy0 - 1.5, f"{gv:g}", size=3.4, color=MUTED, weight="normal", va="bottom")
    ax.plot([px(1.0), px(1.0)], [yy0 - 4, yy0 - 46], color=MUTED, lw=1.6,
            alpha=0.8, zorder=2, ls=(0, (4, 3)))

    for i, (lab, val, note, col, grow) in enumerate(rows):
        v = float(val)
        y = yy0 - 12 - i * 17
        loc = seg(t, 0.4 + i * 0.45, 1.5 + i * 0.45)
        xa, xb = px(1.0), px(v)
        cur = xa + (xb - xa) * loc
        ax.plot([min(xa, cur), max(xa, cur)], [y, y], color=col, lw=5,
                solid_capstyle="round", alpha=0.9, zorder=4)
        ax.add_patch(Circle((cur, y), 2.6, facecolor=col, edgecolor="#0B0E14",
                            lw=1.2, zorder=5))
        txt(ax, x0 - 3, y + 1.8, lab, size=4.0, ha="right", alpha=loc)
        txt(ax, x0 - 3, y - 3.4, note, size=3.0, color=MUTED, weight="normal",
            ha="right", alpha=loc)
        txt(ax, x1 + 2, y, f"{val}×", size=5.2, color=col, ha="left", alpha=loc)

    a = seg(t, 1.9, 2.6)
    txt(ax, 54, 90, "年龄段是最强的风险因子", size=4.8, color=RED, alpha=a)
    txt(ax, 54, 81, "而从墙底摔，重伤风险只有从墙顶摔的 18%", size=4.0, alpha=a)
    txt(ax, 54, 72, "换句话说：离垫子近一点，反而更安全", size=3.6,
        color=MUTED, weight="normal", alpha=a)
    b = seg(t, 2.6, 3.3)
    txt(ax, 54, 58, "数据来源：Beurienne et al. 2025, Front. Sports Act. Living",
        size=2.8, color=MUTED, weight="normal", alpha=b)


def shot6(ax, t, dur):
    """收尾：动作示意 + CTA"""
    grad = np.linspace(0, 1, 256).reshape(-1, 1) / 255.0
    rgb = np.zeros((256, 1, 3))
    rgb[:, 0, 0] = 0.42 + 0.22 * grad[:, 0]
    rgb[:, 0, 1] = 0.50 + 0.28 * grad[:, 0]
    rgb[:, 0, 2] = 0.72 + 0.28 * grad[:, 0]
    ax.imshow(rgb, extent=[0, 108, 0, 192], aspect="auto", zorder=0, alpha=0.22)

    a = seg(t, 0.2, 0.9)
    txt(ax, 54, 180, "记住这两点", size=7.2, alpha=a)

    items = [
        ("双手交叉护住前身", "别用手掌去撑", ORANGE),
        ("落地顺势后仰翻滚", "别用脚硬扛", BLUE),
    ]
    for i, (main, sub, col) in enumerate(items):
        y = 118 - i * 42
        loc = seg(t, 0.7 + i * 0.45, 1.6 + i * 0.45)
        panel(ax, 12, y, 84, 32, fc=PANEL, ec=col, lw=1.4, alpha=0.25 + 0.75 * loc, zorder=3)
        ax.add_patch(Circle((24, y + 16), 6, facecolor=col, edgecolor="none",
                            alpha=0.18 + 0.5 * loc, zorder=4))
        txt(ax, 24, y + 16, str(i + 1), size=4.6, color="#0B0E14", alpha=loc, zorder=6)
        txt(ax, 36, y + 21, main, size=4.8, ha="left", alpha=loc)
        txt(ax, 36, y + 11, sub, size=3.4, color=MUTED, weight="normal",
            ha="left", alpha=loc)

    b = seg(t, dur * 0.62, dur * 0.78)
    panel(ax, 12, 34, 84, 40, fc=PANEL2, ec=ORANGE, lw=1.6, alpha=b, zorder=3)
    txt(ax, 54, 65, "完整解读", size=5.0, color=ORANGE, alpha=b)
    txt(ax, 54, 54, "公众号搜：抱石怎么摔才不伤", size=4.2, alpha=b)
    txt(ax, 54, 45, "245 人 · 301 起损伤 · 法国队列", size=3.4, color=MUTED,
        weight="normal", alpha=b)
    txt(ax, 54, 38, "本文为研究解读，非医疗建议", size=3.0, color=MUTED,
        weight="normal", alpha=seg(t, dur * 0.85, dur * 0.95))


SHOTS = [(shot1_kenburns if S1_VARIANT == "kenburns" else shot1),
         shot2, shot3, shot4, shot5, shot6]


# ---------------------------------------------------------------- timing
def wav_dur(p):
    with wave.open(p, "rb") as w:
        return w.getnframes() / w.getframerate()


def build_timeline():
    durs = [wav_dur(os.path.join(VO_DIR, f"vo_{i}.wav")) for i in range(1, 7)]
    starts, t = [], 0.0
    for d in durs:
        starts.append(t)
        t += d + GAP
    total = starts[-1] + durs[-1] + TAIL
    return durs, starts, total


# ---------------------------------------------------------------- main
def render_preview():
    os.makedirs(OUT, exist_ok=True)
    durs, starts, total = build_timeline()
    for i in range(6):
        fig, ax = new_fig()
        SHOTS[i](ax, durs[i] * 0.5, durs[i])
        cues = cue_times(split_cues(VO_TEXT[i + 1]), durs[i])
        draw_subtitle(ax, cues, durs[i] * 0.5, durs[i])
        p = os.path.join(OUT, f"preview_s{i+1}.png")
        fig.savefig(p, dpi=100, facecolor=BG)
        plt.close(fig)
        print("preview:", p)
    print("durs:", [f"{d:.2f}" for d in durs], "total:", f"{total:.2f}s")


def render_full():
    os.makedirs(OUT, exist_ok=True)
    durs, starts, total = build_timeline()
    print("shot durs:", [f"{d:.2f}" for d in durs])
    print("starts   :", [f"{s:.2f}" for s in starts])
    print("total    :", f"{total:.2f}s")

    cues_all = []
    for i in range(6):
        cues_all += [(c0 + starts[i], c1 + starts[i], c) for c0, c1, c
                     in cue_times(split_cues(VO_TEXT[i + 1]), durs[i])]

    silent = os.path.join(OUT, f"_silent{RENDER_TAG}.mp4")
    final = os.path.join(OUT, f"抱石怎么摔才不伤{RENDER_TAG}.mp4")

    cmd = [
        FFMPEG, "-y", "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", f"{W}x{H}",
        "-r", str(FPS), "-i", "-",
        "-an", "-c:v", "libx264", "-preset", "medium", "-crf", "20",
        "-pix_fmt", "yuv420p", silent,
    ]
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE)

    nframes = int(round(total * FPS))
    for fi in range(nframes):
        t = fi / FPS
        si = 0
        for k in range(6):
            if t >= starts[k]:
                si = k
        local = max(0.0, t - starts[si])
        local = min(local, durs[si])
        fig, ax = new_fig()
        SHOTS[si](ax, local, durs[si])
        draw_subtitle(ax, cues_all, t, total)
        # 直接取 canvas 的 RGBA 缓冲并丢掉 alpha —— savefig(format="raw") 会输出
        # 四通道 RGBA，喂给 -pix_fmt rgb24 会让每帧错位 33%（画面变彩色条纹）
        fig.canvas.draw()
        frame = np.asarray(fig.canvas.buffer_rgba())[:, :, :3]
        plt.close(fig)
        proc.stdin.write(frame.tobytes())
        if fi % 150 == 0:
            print(f"  frame {fi}/{nframes}", flush=True)
    proc.stdin.close()
    proc.wait()
    print("silent video ->", silent)

    # 音轨：按 starts 把各 vo 摆到对应位置
    audio = os.path.join(OUT, f"_voice{RENDER_TAG}.wav")
    sr = 24000
    track = np.zeros(int(total * sr) + sr, dtype=np.int16)
    for i in range(6):
        p = os.path.join(VO_DIR, f"vo_{i+1}.wav")
        with wave.open(p, "rb") as w:
            srate = w.getframerate()
            data = np.frombuffer(w.readframes(w.getnframes()), dtype=np.int16)
        if srate != sr:
            idx = np.linspace(0, len(data) - 1, int(len(data) * sr / srate))
            data = np.interp(idx, np.arange(len(data)), data).astype(np.int16)
        off = int(starts[i] * sr)
        n = min(len(data), len(track) - off)
        track[off:off + n] = data[:n]
    with wave.open(audio, "wb") as w:
        w.setnchannels(1)
        w.setsampwidth(2)
        w.setframerate(sr)
        w.writeframes(track[: int(total * sr)].tobytes())
    print("audio ->", audio)

    mux = [
        FFMPEG, "-y", "-i", silent, "-i", audio,
        "-c:v", "copy", "-c:a", "aac", "-b:a", "160k", "-shortest",
        "-movflags", "+faststart", final,
    ]
    subprocess.run(mux, check=True)
    print("FINAL ->", final)


def render_smoke(nframes=60):
    """只渲前 nframes 帧做管线自检：验证 RGB 通道数正确（抽帧不能是条纹）。"""
    os.makedirs(OUT, exist_ok=True)
    durs, starts, total = build_timeline()
    path = os.path.join(OUT, f"_smoke{RENDER_TAG}.mp4")
    cmd = [FFMPEG, "-y", "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", f"{W}x{H}",
           "-r", str(FPS), "-i", "-", "-an", "-c:v", "libx264", "-preset", "ultrafast",
           "-crf", "18", "-pix_fmt", "yuv420p", path]
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE)
    for fi in range(nframes):
        t = fi / FPS
        fig, ax = new_fig()
        shot1 = SHOTS[0]
        shot1(ax, t, durs[0])
        txt(ax, 54, 22, "管线自检 smoke frame %d" % fi, size=5.6, color="#FFFFFF", stroke=5)
        fig.canvas.draw()
        frame = np.asarray(fig.canvas.buffer_rgba())[:, :, :3]
        plt.close(fig)
        proc.stdin.write(frame.tobytes())
    proc.stdin.close()
    proc.wait()
    print("smoke ->", path, f"({nframes} frames, {nframes/FPS:.1f}s)")
    chk = os.path.join(OUT, f"_smoke_check{RENDER_TAG}.png")
    subprocess.run([FFMPEG, "-v", "error", "-y", "-i", path,
                    "-frames:v", "1", chk], check=True)
    print("check ->", chk)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "preview"
    if mode == "preview":
        render_preview()
    elif mode == "smoke":
        render_smoke()
    else:
        render_full()
