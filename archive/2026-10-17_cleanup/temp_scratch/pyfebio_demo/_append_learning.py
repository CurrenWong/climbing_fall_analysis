"""把 pyfebio/FEBio 相关的发现追加到 .learnings/ERRORS.md（逐条独立、幂等）。"""

from pathlib import Path

P = Path(__file__).resolve().parents[2] / ".learnings" / "ERRORS.md"

ADD_006 = r"""
---

## [ERR-20261003-006] pyfebio 0.3.0 能写不能读：`Literal[0,1]` 字段反序列化必炸

**类型**：第三方库 bug（往返失败）
**发现方式**：官方 README 样例跑完，用 `Model.from_xml()` 读回自己刚写的 .feb
**状态**：未修（上游 bug）—— 用法上按「只写不读」处理

### 症状

`Model.from_xml(feb_bytes)` 抛 77 个 ValidationError：
`boundary_.all_bcs.0.BCZeroDisplacement.x_dof: Input should be 0 or 1, input_type=str`

### 根因（最小复现见 `temp/pyfebio_demo/repro.py`）

pydantic 的 `Literal[0, 1]` **只接受 Python int，不接受 str `"1"`**：

| 输入 | `Literal[0,1]` | `int` |
|---|---|---|
| `1` (int) | OK | OK |
| `"1"` (str) | **FAIL** | OK |

而 pydantic-xml 反序列化时给的**永远是 str**。所以
`pyfebio/boundary.py` 里所有 `x_dof/y_dof/z_dof/relative/aggressiveness/
dtforce/cycle_buffer` 这些 `Literal[0,1]` 字段读回来必然失败。
`control_.solver.*.qn_method` 的 `type` 判别失败会连锁炸 20+ 个字段。

### 附带的调用坑

- `Model.from_xml(str)` → `ValueError: Unicode strings with encoding declaration
  are not supported`（.feb 带 `<?xml encoding='ISO-8859-1'?>`）→ **必须传 bytes**
- `Model.parse_file(path)` 也不行 —— 那是 pydantic 的 **JSON** 解析，直接
  `JSONDecodeError`（且已 deprecated）

### 正确读法

```python
m = pyfebio.model.Model.from_xml(Path("x.feb").read_bytes())  # bytes
```

但即使传了 bytes 仍会撞上面的 Literal bug。

### 结论

pyfebio 0.3.0 是 **generator-only** 的：写 .feb 可用，读回来不可用。
依赖它的流程只能是「Python 对象 → .feb → febio4」单向。
另外 `run_model()` 只是 `subprocess.run("febio4 ...")`，不在 PATH 时
返回 rc=1 且只打印一行 cmd 错误，**不抛异常**。

**关联**：`temp/pyfebio_demo/{demo,roundtrip,repro}.py`
"""

ADD_007 = r"""
---

## [ERR-20261003-007] FEBio 4.13.0 不在系统 PATH；XPLT→HDF5 的 dataset 布局反直觉

**类型**：环境 + 集成坑
**状态**：已绕过（样例已跑通并通过物理校验）

### 坑 1：febio4 装在 D 盘但没进 PATH

安装位置 `D:\Program\FEBioStudio\bin\febio4.exe`（FEBio 4.13.0），
**不在系统 PATH 里**，也不在用户 PATH 里。每次跑都得手动：

```cmd
set "PATH=D:\Program\FEBioStudio\bin;%PATH%" && python demo.py --run
```

`run_model()` 是裸 `subprocess.run("febio4 -i ...")`，不做路径探测，
所以没 prepend 就静默返回 rc=1（**只打一行 cmd 错误，不抛异常**）。

### 坑 2：xplt.to_hdf5() 的 dataset 布局

`pyfebio.xplt` 需要**显式 import**（`__init__.py` 没导出），
转换入口是 **`xplt.to_hdf5(src, dst)`**（不是 `convert_xplt_to_hdf5`）。

转换后 layout（踩了两次才摸清）：

| 路径 | dtype | shape | 备注 |
|---|---|---|---|
| `meshes/0/nodes` | 结构化 `[id,x,y,z]` | (8,) | **是** recarray |
| `states/N/node_data/displacement/nodes` | `float32` | (8,3) | **不是**结构化，**列=xyz 分量** |
| `states/N/element_data/stress/box` | `float32` | (1,6) | Voigt 序 xx yy zz xy yz zx |

两个关键点：

1. **每个 state 是独立 dataset**（`states/0/...`、`states/1/...`），
   **没有**堆叠成时间轴的三维数组 → 不能 `f[key][:]` 一次拿全部帧。
2. 位移 dataset 是**纯 float32 (n,3)**，不是 recarray。
   写 `u[4]["z"]` 或 `u["z"][4]` 都会 `IndexError`
   （看着像结构化数组，习惯性写法必踩），正确写法 `u[4, 2]` / `u[:, 2]`。

### 顺带：`cmd.exe` 吞多行 Python

`python -c "..."` 里写多行 + 缩进，cmd.exe 会把换行折成空格导致
`IndentationError: unexpected indent`。**多行 Python 一律写成 `.py` 再跑。**

**关联**：`temp/pyfebio_demo/{demo,read_xplt,probe_h5}.py`
"""

text = P.read_text(encoding="utf-8")
for title, block in (
    ("## [ERR-20261003-006] pyfebio", ADD_006),
    ("## [ERR-20261003-007] FEBio 4.13.0", ADD_007),
):
    if title in text:
        print(f"skip (exists): {title}")
    else:
        P.write_text(P.read_text(encoding="utf-8").rstrip() + "\n" + block, encoding="utf-8")
        print(f"appended: {title}  ({P.stat().st_size} bytes total)")
