from pathlib import Path

M = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md")
L = 10240
s = M.read_text(encoding="utf-8")
CUTS = [
    ("| WXML 报 `unexpected token .` | WXML 表达式**不能调方法、不能用 `?.`** ⇒ 格式化全在 JS 层做成 `xxxText` |",
     "| WXML `unexpected token .` | WXML 表达式**不能调方法/不能 `?.`** ⇒ 格式化在 JS 层做 `xxxText` |"),
    ("| 小程序 input 高度比 picker 矮 | 原生组件不响应 `line-height` ⇒ 用 `min-height: 88rpx` |",
     "| 小程序 input 比 picker 矮 | 原生组件不响应 `line-height` ⇒ `min-height: 88rpx` |"),
    ("| miniprogram-ci `Unexpected token .` | babel 不支持 `?.` ⇒ `(a && a.b)` 或开 swc |",
     "| miniprogram-ci `Unexpected token .` | babel 不支持 `?.` ⇒ `(a && a.b)` |"),
    ("| 金额 `151.2699999998` | JS 浮点累加 ⇒ 全局 `round2`（`Number.EPSILON` 不能省）|",
     "| 金额 `151.2699999998` | JS 浮点累加 ⇒ 全局 `round2`（`EPSILON` 不能省）|"),
]
for a, b in CUTS:
    if a in s:
        s = s.replace(a, b, 1)
    M.write_text(s, encoding="utf-8")
sz = len(s.encode())
print(f"disk = {sz:,} 余量 {L - sz:+,}")
