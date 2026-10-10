"""B1 —— 论文 v2 跨位置数字一致性核对。

两部分：
  P1 精选关键量：每量一个 regex + 期望值，列出所有出现位置与取值，
     多值即 FAIL（除非在 EXPECTED 的白名单内）。
  P2 关键词邻域扫描：对若干关键词，取其 ±N 字符内的数字，
     若同一关键词在不同章节邻域数字集合不同 => 列为「待人工复核」。
P2 只报候选不下结论（数字不同未必错，可能本就不同量）。
"""
from __future__ import annotations

import pathlib
import re
from collections import defaultdict

ROOT = pathlib.Path(__file__).resolve().parents[2]
PAPER = ROOT / "docs/paper/论文v2_数学建模版.md"

# --------------------------------------------------------------------------
# 分节
# --------------------------------------------------------------------------
HEAD = re.compile(r"^(#{1,3})\s+(.+)$", re.M)


def split_sections(text: str) -> list[tuple[str, str]]:
    marks = [(m.start(), m.group(2).strip()) for m in HEAD.finditer(text)]
    out = []
    for i, (pos, title) in enumerate(marks):
        end = marks[i + 1][0] if i + 1 < len(marks) else len(text)
        out.append((title, text[pos:end]))
    return out


# --------------------------------------------------------------------------
# P1 精选关键量
# --------------------------------------------------------------------------
# name: (regex, 期望值白名单 —— 出现的每个数值都必须在内；None=只报不判)
CURATED = {
    "踝损伤率 34.5%": (r"34\.5\s*\\?%", {"34.5"}),
    "踝骨折率 22.8%": (r"22\.8\s*\\?%", {"22.8"}),
    "临床踝损伤 40%": (r"\b40(\.0)?\s*\\?%", {"40"}),
    "临床踝扭伤 28%": (r"\b28(\.0)?\s*\\?%", {"28"}),
    "临床踝骨折 8%": (r"对\s*8\s*\\?%", {"8"}),
    "倍数 1.16x": (r"1\.16\s*[x×]", {"1.16"}),
    "倍数 2.85x": (r"2\.85\s*[x×]", {"2.85"}),
    "撕脱阈 215 N": (r"\b215\b", {"215"}),
    "撕脱旧阈 300": (r"\b300(\.0)?\s*N", {"300"}),
    "alpha=0.477": (r"0\.477", {"0.477"}),
    "alpha_post=0.4800": (r"0\.4800", {"0.48"}),
    "V3 阈 12.5 ms": (r"12\.5\s*ms", {"12.5"}),
    "V3 阈 0.0125 s": (r"0\.0125", {"0.0125"}),
    "V3【垫】下界 26.6": (r"\b26\.6\b", {"26.6"}),
    "V3【垫】上界 37.0": (r"\b37\.0\b", {"37.0"}),
    "V8 阈 5 ms@20m": (r"≤\s*5\s*ms|5\s*ms\s*@|0\.005", None),
    "V8b【刚】14.0": (r"\b14\.0\b", {"14.0"}),
    "V8a 起 18.8": (r"\b18\.8\b", {"18.8"}),
    "S5 34.4 kN": (r"34\.4", {"34.4"}),
    "S5 0.574 MPa": (r"0\.574", {"0.574"}),
    "撤销 2.70x": (r"2\.70\s*[x×]", {"2.70"}),
    "撤销 15.14 cm": (r"15\.14", {"15.14"}),
    "k_neck 标定 2500": (r"\b2500\b", {"2500"}),
    "P5 h5@10 11.40": (r"11\.40", {"11.40"}),
    "P5 h2@30 10.00": (r"10\.00", {"10.00"}),
    "P5 h3@30 9.60": (r"9\.60", {"9.60"}),
    "P5 20m a100 5.2": (r"\b5\.2\b", {"5.2"}),
    "P5 20m a300 2.8": (r"\b2\.8\b", {"2.8"}),
    "两因子 单改工况 22.2": (r"\b22\.2\b", {"22.2"}),
    "两因子 单改k 31.4": (r"\b31\.4\b", {"31.4"}),
    "两因子 双改 10.8": (r"\b10\.8\b", {"10.8"}),
    "工况占比 67.6": (r"67\.6", {"67.6"}),
    "利用率位移 -25%": (r"−?\s*25\s*%", {"25"}),
    "util 上界 +27%": (r"\+\s*27\s*%", {"27"}),
    "Q1 峰值力 -24%": (r"−?\s*24\s*\\?%", {"24"}),
    "Q1 接触时长 +139%": (r"139\s*\\?%", {"139"}),
    "无旋转 0/321": (r"0/321", {"0/321"}),
}

# 关键词邻域（P2）
KEYWORDS = ["踝骨折", "踝损伤", "峰值时刻", "冲击窗", "撕脱", "峰值力",
            "利用率", "载荷分享", "k_neck", "距下"]


def main() -> int:
    text = PAPER.read_text(encoding="utf-8")
    secs = split_sections(text)
    print("=" * 84)
    print("P1  精选关键量：跨章节取值是否一致")
    print("=" * 84)
    fails = []
    print(f"{'量':<26} {'出现':>4}  章节分布 / 取值")
    print("-" * 84)
    for name, (pat, allowed) in CURATED.items():
        hits = []
        for title, body in secs:
            for m in re.finditer(pat, body):
                raw = m.group(0).replace("\n", " ")
                hits.append((title, raw))
        if not hits:
            print(f"{name:<26} {'0':>4}  -- 未出现（论文可能不引该量）")
            continue
        vals = {re.sub(r"[^\d./]", "", r) for _, r in hits}
        secs_hit = sorted({t for t, _ in hits})
        flag = ""
        if allowed is not None:
            # 数值归一化比较（'40.0' 与 '40'、'0.4800' 与 '0.48' 同值）
            def to_num(v: str):
                v = v.rstrip(".")
                if not v or re.fullmatch(r"\d+/\d+", v):
                    return v
                try:
                    return round(float(v), 6)
                except ValueError:
                    return v
            norm = {to_num(v) for v in vals if v != ""}
            okset = {to_num(a) for a in allowed}
            bad = {v for v in norm if v not in okset}
            if bad:
                flag = "  <<< MISMATCH " + repr(bad)
                fails.append((name, bad))
        print(f"{name:<26} {len(hits):>4}  {','.join(secs_hit[:4])}"
              f"{'…' if len(secs_hit) > 4 else ''}  vals={sorted(vals)}{flag}")

    print()
    print("=" * 84)
    print("P2  同关键词跨章节「接近但不相等」的数字 = 舍入/打架候选")
    print("=" * 84)
    # 过滤：年份、章节号（如 4.6）、一位数
    def keep(v: float) -> bool:
        if v < 10:
            return False
        if 1900 <= v <= 2100:
            return False
        return True

    cands = []
    ctx_of: dict = {}
    for kw in KEYWORDS:
        locs = defaultdict(dict)          # section -> {value: snippet}
        for title, body in secs:
            if kw not in body:
                continue
            for m in re.finditer(re.escape(kw), body):
                ctx = body[max(0, m.start() - 45):m.end() + 45].replace("\n", " ")
                for n in re.findall(r"\d+\.\d+|\d+", ctx):
                    try:
                        f = float(n)
                    except ValueError:
                        continue
                    if keep(f):
                        locs[title].setdefault(f, ctx[:96])
                        ctx_of[(kw, title, f)] = ctx[:96]
        keys = sorted(locs)
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                for a in locs[keys[i]]:
                    for b in locs[keys[j]]:
                        if a == b:
                            continue
                        lo, hi = min(a, b), max(a, b)
                        if hi / lo <= 1.15:      # 接近但不等
                            cands.append((kw, keys[i], a, keys[j], b))
    # 自动判掉明显假警报：
    #  (a) 上下文里该数落在 [NN] 文献引用里
    #  (b) 两侧数字紧邻的**单位不同**（kN / ms / MPa / m / %）=> 本就是不同量
    UNIT = re.compile(r"(kN|MPa|ms|N·m|N/mm|cm|mm|%|m)(?![a-zA-Z])")
    cands_real, auto_dismissed = [], []
    for kw, s1, a, s2, b in cands:
        c1, c2 = ctx_of.get((kw, s1, a), ""), ctx_of.get((kw, s2, b), "")
        # 引用号：在**整节正文**里查（上下文 96 字窗口常把 [22,26,27] 截断）
        body1 = dict(secs).get(s1, "")
        body2 = dict(secs).get(s2, "")
        if float(a).is_integer() and re.search(
                rf"\[[^\]]*\b{int(a)}\b[^\]]*\]", body1):
            auto_dismissed.append((kw, s1, a, s2, b, "citation"))
            continue
        if float(b).is_integer() and re.search(
                rf"\[[^\]]*\b{int(b)}\b[^\]]*\]", body2):
            auto_dismissed.append((kw, s1, a, s2, b, "citation"))
            continue
        # 年份片段（2026-10-08 -> 202 / 08）
        if (1900 <= a <= 2100 and "2026" in c1.replace("202", "2026")):
            auto_dismissed.append((kw, s1, a, s2, b, "date"))
            continue
        u1 = set(UNIT.findall(c1))
        u2 = set(UNIT.findall(c2))
        if u1 and u2 and not (u1 & u2):
            auto_dismissed.append((kw, s1, a, s2, b,
                                   f"units {'/'.join(u1)} vs {'/'.join(u2)}"))
            continue
        cands_real.append((kw, s1, a, s2, b))

    seen, shown = set(), 0
    for kw, s1, a, s2, b in cands_real:
        key = (kw, round(min(a, b), 3), round(max(a, b), 3))
        if key in seen:
            continue
        seen.add(key)
        shown += 1
        print(f"  [{kw}] {s1}={a}  vs  {s2}={b}   (差 {abs(a-b):.3g})")
        print(f"        A: ...{ctx_of.get((kw, s1, a), '')}...")
        print(f"        B: ...{ctx_of.get((kw, s2, b), '')}...")
    print()
    print(f"  自动判掉（假警报）: {len(auto_dismissed)}")
    rs = defaultdict(int)
    for *_, why in auto_dismissed:
        rs[why if not why.startswith(("units", "citation")) else
           why.split()[0] if why.startswith("units") else "citation"] += 1
    for k, v in sorted(rs.items(), key=lambda x: -x[1]):
        print(f"     {k:<16} {v}")
    print(f"  **需人工复核**: {len(seen)} 组"
          f"  （以上列全部；上方未出现的即已自动判掉）")
    print("=" * 84)
    print("RESULT  " + ("PASS" if not fails
                        else "FAIL: " + repr([f[0] for f in fails])))
    print("=" * 84)
    return 0 if not fails else 1


if __name__ == "__main__":
    raise SystemExit(main())
