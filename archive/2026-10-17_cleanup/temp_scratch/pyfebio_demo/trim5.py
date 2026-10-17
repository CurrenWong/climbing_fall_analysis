"""MEMORY 精修：削 371+ B。"""
import io

P = r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md"
t = io.open(P, encoding="utf-8").read()
print(f"改前 {len(t.encode('utf-8'))} B")
REPL = [
    ("- 🧭 **`dispx` 下 `DISP_MM>0` = 压实**（`[ERR-20261006-003]` 三证：`.feb` 写 z=−DISP｜"
     "dz 朝距骨｜接触 149.6 vs 0 N）⇒ `+1.6` 的 **592.6 N 是真压实**；"
     "`[ERR-20261005-009]` 作废；`LOAD_SIGN` 对 `dispx` 无效",
     "- 🧭 **`dispx` 下 `DISP_MM>0` = 压实**（`[ERR-20261006-003]`；`.feb` z=−DISP｜dz 朝距骨｜"
     "接触 149.6 vs 0 N）⇒ `+1.6` 的 **592.6 N 是真压实**；`[ERR-20261005-009]` 作废"),
    ("- ✅ **压缩曲线**：+0.25/0.5/1.0/**1.6**/2.5 mm → 124.0/230.7/418.1/**592.6**/615.4 N"
     "⇒ 600 N @≈+1.9mm｜`[LRN-072]` 负值序列是拉伸、作废",
     "- ✅ **压缩曲线**：+0.25/0.5/1.0/**1.6**/2.5mm → 124/230.7/418.1/**592.6**/615.4 N"
     "⇒ 600N @≈+1.9mm｜`[LRN-072]` 负值=拉伸、作废"),
    ("- 🎯 **骨折裕度（P1 首版）**：@592.6N 腓骨皮质 **11.40 MPa/裕度 0.38**（≈断 1560 N）"
     "> 距骨 6.57 > 跟骨 6.21 > 胫骨 5.97，全 ≪1；**腓骨最受力**（合临床）",
     "- 🎯 **骨折裕度（P1 首版）**：@592.6N 腓骨皮质 **11.40 MPa/裕度 0.38**（≈断 1560N）"
     "> 距骨 6.57 > 跟骨 6.21 > 胫骨 5.97，全 ≪1；**腓骨最受力**"),
    ("- ✅ **P1 应力输入已独立验证**（位移→梯度→本构）：**tet4 域比值 0.98~1.03** ⇒ h5 骨应力可信；"
     "`navicular_spon` 比值 **9.6** ⇒ **该域勿用**；σ_max 104.855 = 伪影（重构 8.07）",
     "- ✅ **P1 应力输入已独立验证**（位移→梯度→本构）：**tet4 比值 0.98~1.03** ⇒ h5 骨应力可信；"
     "`navicular_spon` 比值 **9.6** ⇒ **勿用**；σ_max 104.855=伪影（重构 8.07）"),
    ("- 🔬 **软骨网格**：① `area_proj` **quad4 翘曲 bug 已修**（47.8→**555.5** mm²）"
     "② **11 个翻转单元全在关节面⇒不可剔除**；`detJ_penta6_full` **只诊断**（`[LRN-087]`）"
     "③ `cart_m1_fix` 有**贴合度回归**（+5° 张开 0.598mm）⇒ **保持 `cart_m1` 为基线**",
     "- 🔬 **软骨网格**：① `area_proj` **quad4 翘曲 bug 已修**（47.8→**555.5**mm²）"
     "② **11 个翻转单元全在关节面⇒不可剔除**；`detJ_penta6_full` **只诊断**（`[LRN-087]`）"
     "③ `cart_m1_fix` **贴合度回归**（+5° 张开 0.598mm）⇒ **保持 `cart_m1` 基线**"),
    ("- 🧨 **坑**：`lc`→`<scale>`｜`sym_stiffness` 在 `<solver>`｜删旧 `.feb`+查 rc｜"
     "**patch test 须用畸变单元**｜**h5 nodesets ≠ `.feb` NodeSet**（差 46/35/80/59）",
     "- 🧨 **坑**：`lc`→`<scale>`｜`sym_stiffness` 在 `<solver>`｜删旧 `.feb`+查 rc｜"
     "**patch test 须用畸变单元**｜**h5 nodesets ≠ `.feb` NodeSet**（差 46/35/80/59）"),
]
hit = 0
for a, b in REPL:
    if a in t:
        t = t.replace(a, b); hit += 1
    else:
        print("  ⚠ 未命中：", a[:34].replace("\n", "\\n"))
# 兜底：若仍超，削 Anderson 行的冗余
if len(t.encode("utf-8")) > 10240:
    A = ("- **Anderson 真值**（`[LRN-067]`）：all-hex｜距骨**弱弹簧 500 N/mm**+解除内外翻内旋｜"
         "**骨刚性+reference node**｜不建韧带｜E=12 ν=0.42｜落位 200 N｜精度 10~15%")
    B = ("- **Anderson 真值**（`[LRN-067]`）：all-hex｜距骨**弱弹簧 500N/mm**+解除内外翻内旋｜"
         "**骨刚性+reference node**｜不建韧带｜E=12 ν=0.42｜落位 200N")
    if A in t:
        t = t.replace(A, B)
io.open(P, "w", encoding="utf-8").write(t)
n = len(t.encode("utf-8"))
print(f"命中 {hit}/{len(REPL)}；改后 {n} B  "
      + ("✅ 余量 %d" % (10240-n) if n <= 10240 else "❌ 超 %d" % (n-10240)))
