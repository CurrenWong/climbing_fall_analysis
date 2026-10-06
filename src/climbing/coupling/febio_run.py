"""FEBio 求解器定位与安全调用（方案 S0.2）。

背景
----
``pyfebio.model.run_model()`` 本质是裸 ``subprocess.run("febio4 ...")``：
FEBio 装在 ``D:\\Program\\FEBioStudio\\bin`` 但**不在 PATH** 时，它
**静默返回 rc=1**，只打一行 cmd 错误、不抛异常
（见 ``.learnings/ERRORS.md`` 的 ``[ERR-...-007]``）。

本模块把"找 exe → 跑 → 检查退出码/输出"显式化，使这类静默失效**不可能通过**：
找不到就抛 :class:`FebioNotFound`，跑挂就抛 :class:`FebioRunError`（带日志尾部）。

用法
----
.. code-block:: python

    from climbing.coupling.febio_run import find_febio, run_febio

    exe = find_febio()                       # -> Path；找不到抛 FebioNotFound
    print(probe_febio(exe))                  # "4.13.0"（无参 banner 探测）
    run_febio("bone.feb")                    # rc!=0 / 致命错误 -> 抛 FebioRunError

定位顺序（先命中先用）::

    1. explicit 参数
    2. 环境变量 FEBIO_EXE
    3. 已知安装位置 DEFAULT_CANDIDATES
    4. PATH 上的 febio4 / febio4.exe / febio
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

#: 已知安装位置（按优先级）。可用环境变量 ``FEBIO_EXE`` 覆盖。
DEFAULT_CANDIDATES: tuple[str, ...] = (
    r"D:\Program\FEBioStudio\bin\febio4.exe",
)

#: PATH 上可接受的可执行名
PATH_NAMES: tuple[str, ...] = ("febio4", "febio4.exe", "febio")

#: 覆盖用环境变量
ENV_VAR = "FEBIO_EXE"

#: 从 banner 抓版本号（"version 4.13.0"）
_VERSION_RE = re.compile(r"version\s+([0-9][0-9.]*)", re.IGNORECASE)


class FebioNotFound(FileNotFoundError):
    """找不到可运行的 FEBio 可执行文件。"""


class FebioRunError(RuntimeError):
    """FEBio 运行失败（非 0 退出码，或日志里有 FATAL ERROR）。"""


# --------------------------------------------------------------------------
# 探测
# --------------------------------------------------------------------------
def _banner(exe: Path, timeout: float = 30.0) -> str | None:
    """无参调用 exe，返回 stdout+stderr；不可执行时返回 None。

    FEBio 无参时会打印 banner（含版本号）并以非 0 退出，所以这里**不**看退出码，
    只看能否启动且输出像 FEBio。
    """
    try:
        proc = subprocess.run(
            [str(exe)], capture_output=True, text=True, timeout=timeout
        )
    except (OSError, subprocess.SubprocessError, ValueError):
        return None
    return f"{proc.stdout or ''}\n{proc.stderr or ''}"


def probe_febio(exe: str | os.PathLike, *, timeout: float = 30.0) -> str | None:
    """探测 exe 是否为可运行的 FEBio，返回版本号字符串（如 ``"4.13.0"``），否则 None。"""
    out = _banner(Path(exe), timeout=timeout)
    if out is None or "FEBio" not in out:
        return None
    m = _VERSION_RE.search(out)
    return m.group(1) if m else "unknown"


def _not_found_msg(tried: list[tuple[str, Path]], *, explicit: bool) -> str:
    lines = [
        "找不到可运行的 FEBio 可执行文件。" if not explicit else "指定的 FEBio 路径不可用。",
        "",
        "已尝试：",
    ]
    lines += [f"  - [{src}] {p}" for src, p in tried] or ["  （无）"]
    lines += [
        "",
        "修复方式（任选其一）：",
        f"  - 设环境变量 {ENV_VAR} 指向 febio4.exe；",
        "  - 把 FEBio 的 bin 目录加入 PATH；",
        "  - 调用时显式传入 exe 路径。",
        "",
        "注意：pyfebio 的 run_model() 在找不到 FEBio 时只会静默返回 rc=1，",
        "本函数的存在就是为了让这种情况抛错而不是被忽略。",
    ]
    return "\n".join(lines)


def find_febio(
    explicit: str | os.PathLike | None = None,
    *,
    verify: bool = False,
) -> Path:
    """定位 FEBio 可执行文件，返回绝对 :class:`Path`。

    Parameters
    ----------
    explicit:
        直接指定的路径；给了但不可用则立即抛错（不会回退到别的候选）。
    verify:
        为 True 时对候选做一次无参 banner 探测，只有**真的能启动**才算命中。

    Raises
    ------
    FebioNotFound
        所有候选都不可用。
    """
    tried: list[tuple[str, Path]] = []

    def _accept(path: Path, source: str) -> Path | None:
        tried.append((source, path))
        if not path.is_file():
            return None
        if verify and probe_febio(path) is None:
            return None
        return path.resolve()

    # 1) 显式传入
    if explicit is not None:
        p = _accept(Path(explicit), "explicit")
        if p is None:
            raise FebioNotFound(_not_found_msg(tried, explicit=True))
        return p

    # 2) 环境变量
    env = os.environ.get(ENV_VAR)
    if env:
        p = _accept(Path(env), f"env:{ENV_VAR}")
        if p is not None:
            return p

    # 3) 已知安装位置
    for cand in DEFAULT_CANDIDATES:
        p = _accept(Path(cand), "candidate")
        if p is not None:
            return p

    # 4) PATH
    for name in PATH_NAMES:
        found = shutil.which(name)
        if found:
            p = _accept(Path(found), f"PATH:{name}")
            if p is not None:
                return p

    raise FebioNotFound(_not_found_msg(tried, explicit=False))


# --------------------------------------------------------------------------
# 运行
# --------------------------------------------------------------------------
def _tail(text: str, n: int = 40) -> str:
    lines = text.splitlines()
    return "\n".join(lines[-n:]) if len(lines) > n else text


def _log_tail(feb_path: Path, cwd: Path, n: int = 40) -> str:
    """FEBio 的 .log 通常写在输入文件旁；也兜底看 cwd。"""
    for cand in (feb_path.with_suffix(".log"), cwd / (feb_path.stem + ".log")):
        if cand.is_file():
            try:
                return _tail(cand.read_text(errors="replace"), n)
            except OSError:
                pass
    return ""


def run_febio(
    feb_path: str | os.PathLike,
    *,
    workdir: str | os.PathLike | None = None,
    exe: str | os.PathLike | None = None,
    silent: bool = False,
    timeout: float | None = None,
    check_fatal: bool = True,
) -> int:
    """运行 FEBio 求解一个 ``.feb``，返回退出码（成功即 0）。

    与 ``pyfebio.model.run_model()`` 的区别：**失败会抛错，绝不静默返回非 0**。

    Parameters
    ----------
    feb_path:
        输入 ``.feb`` 文件。
    workdir:
        工作目录（FEBio 的输出相对此目录）；默认 = 输入文件所在目录。
    exe:
        显式指定 FEBio 可执行文件；默认走 :func:`find_febio`。
    silent:
        是否加 ``-silent``（抑制控制台输出；错误信息仍会从 .log 兜底）。
    timeout:
        ``subprocess`` 超时（秒）。
    check_fatal:
        退出码为 0 时仍扫描输出里的 ``FATAL ERROR``（防御 FEBio 的罕见情形）。

    Raises
    ------
    FebioNotFound, FileNotFoundError, FebioRunError
    """
    exe_p = find_febio(exe)
    feb_p = Path(feb_path)
    if not feb_p.is_file():
        raise FileNotFoundError(f"输入文件不存在：{feb_p}")

    cwd = Path(workdir) if workdir is not None else feb_p.parent
    cmd = [str(exe_p), "-i", str(feb_p)]
    if silent:
        cmd.append("-silent")

    proc = subprocess.run(
        cmd, capture_output=True, text=True, cwd=str(cwd), timeout=timeout
    )
    combined = f"{proc.stdout or ''}\n{proc.stderr or ''}"

    def _fail(reason: str) -> FebioRunError:
        msg = [
            reason,
            f"命令: {' '.join(cmd)}",
            f"工作目录: {cwd}",
            f"退出码: {proc.returncode}",
        ]
        tail = _log_tail(feb_p, cwd) or _tail(combined)
        if tail.strip():
            msg += ["", "---- 输出（尾部）----", tail, "---------------------"]
        return FebioRunError("\n".join(msg))

    if proc.returncode != 0:
        raise _fail("FEBio 运行失败（非 0 退出码）。")
    if check_fatal and "FATAL ERROR" in combined.upper():
        raise _fail("FEBio 退出码为 0，但输出含 FATAL ERROR。")

    return proc.returncode


# --------------------------------------------------------------------------
# 自检：python -m climbing.coupling.febio_run
# --------------------------------------------------------------------------
def _main() -> int:
    try:
        exe = find_febio(verify=True)
    except FebioNotFound as exc:
        print(f"[FAIL] {exc}")
        return 1
    ver = probe_febio(exe)
    print(f"[ok] FEBio 可执行文件: {exe}")
    print(f"[ok] 版本探测: {ver}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
