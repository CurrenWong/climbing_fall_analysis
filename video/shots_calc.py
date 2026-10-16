# -*- coding: utf-8 -*-
"""按实测语速重算逐镜旁白时长。音节规则：汉字 1 音节，数字按中文读法展开。"""

DIGIT = {"0": "零", "1": "一", "2": "二", "3": "三", "4": "四",
         "5": "五", "6": "六", "7": "七", "8": "八", "9": "九"}

def num_cn(s):
    """阿拉伯数字串 -> 中文读法（够用的近似：整数/一位小数/百分号）。"""
    if isinstance(s, int):
        n = s
    else:
        s = s.strip()
        if "." in s:
            a, b = s.split(".", 1)
            return (num_cn(a) + "点" + "".join(DIGIT[c] for c in b)) if a else ("零点" + "".join(DIGIT[c] for c in b))
        n = int(s)
    # 处理 0-9999
    units = ["", "十", "百", "千"]
    if n == 0:
        return "零"
    if n < 10:
        return DIGIT[str(n)]
    if n < 100:
        t, r = divmod(n, 10)
        head = "十" if t == 1 else DIGIT[str(t)] + "十"
        return head + (DIGIT[str(r)] if r else "")
    if n < 1000:
        t, r = divmod(n, 100)
        return DIGIT[str(t)] + "百" + ("零" + DIGIT[str(r)] if 0 < r < 10 else (num_cn(r) if r else ""))
    t, r = divmod(n, 1000)
    return DIGIT[str(t)] + "千" + ("零" + num_cn(r) if 0 < r < 100 else (num_cn(r) if r else ""))

def syllables(t):
    total, i, marks = 0, 0, 0
    while i < len(t):
        c = t[i]
        if c.isdigit():
            j = i
            while j < len(t) and (t[j].isdigit() or (t[j] == "." and j + 1 < len(t) and t[j + 1].isdigit())):
                j += 1
            chunk = t[i:j]
            if j < len(t) and t[j] == "%":
                chunk += "%"
                j += 1
            spoken = num_cn(chunk) if not chunk.endswith("%") else "百分之" + num_cn(chunk[:-1])
            total += len(spoken)
            i = j
            continue
        if "一" <= c <= "鿿":
            total += 1
        elif c in "，、。：；？！（）":
            marks += 1
        i += 1
    return total, marks

VO = {
    "S1": "抱石的人，迟早要摔一次。可你真的知道自己是怎么摔的吗？",
    "S2": "一项法国研究追踪了245名抱石者的301起损伤：67%伤在下肢，踝关节占四成，伤后不能攀爬的中位时间是90天。",
    "S3": "更关键的是，坠落的方式分两条完全不同的路。第一条：身体保持垂直、坠落中旋转、最后脚先着地——占大多数，伤在脚踝。脚和踝陷进垫子，踝被迫过度外旋，韧带首当其冲。第二条：起手就是后仰，坠落中大幅转体，落地不稳时本能用手撑垫子——伤在肘和肩。",
    "S4": "最常见的前三种场景里，两种都是垂直落体、脚先着地，加起来超过三成，结果都是下肢伤。",
    "S5": "什么让人伤得更重？28到31岁的重伤风险，是32岁以上人群的3.35倍。从墙底摔，比从墙顶摔低82%。",
    "S6": "所以：双手交叉护住身前，落地顺势后仰翻滚。完整解读，看公众号。",
}

RATES = {"A(冷静客观科普·偏快)": 5.12, "B(干脆利落快评·很快)": 5.53, "C(冷静客观科普·中等偏快)": 4.31}

print("逐镜音节数（数字按中文读法展开）")
print("-" * 74)
tot = 0
for k in ["S1", "S2", "S3", "S4", "S5", "S6"]:
    syl, marks = syllables(VO[k])
    tot += syl
    cells = "  ".join(f"r={r:<4} →{syl/r:5.1f}s" for r in RATES.values())
    print(f"{k}  音节={syl:3d}  标点={marks:2d}   {cells}")
print("-" * 74)
print(f"合计音节 = {tot}")
for name, r in RATES.items():
    print(f"  {name:26s} r={r:4.2f} → 全文 {tot/r:5.1f}s")
