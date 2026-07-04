"""
底层 IO 工具层：统一文件句柄管理、格式探测。

- open_sequence_file(): 上下文管理器，透明支持普通文本与 gzip 压缩文件。
- detect_format(): 根据文件内容（首字符）自动判断 FASTA / FASTQ。
- FileFormat: 文件格式枚举。
"""

from __future__ import annotations

import gzip
import io
from contextlib import contextmanager
from enum import Enum
from pathlib import Path
from typing import IO, Iterator, Union

# 类型别名：可接受字符串路径或 Path 对象
PathLike = Union[str, Path]


class FileFormat(Enum):
    """支持的序列文件格式。"""
    FASTA = "fasta"
    FASTQ = "fastq"
    UNKNOWN = "unknown"


@contextmanager
def open_sequence_file(filepath: PathLike, mode: str = "rt") -> Iterator[IO[str]]:
    """
    上下文管理器：自动判断文件是否为 gzip 压缩，返回文本模式的文件句柄。

    Parameters
    ----------
    filepath : str or Path
        文件路径，支持 .fasta / .fastq / .fa / .fq 及其 .gz 压缩版本。
    mode : str
        打开模式，默认 "rt"（文本读取）。

    Yields
    ------
    IO[str]
        文本模式的文件句柄，可逐行迭代。

    Examples
    --------
    >>> with open_sequence_file("example.fasta.gz") as fh:
    ...     for line in fh:
    ...         print(line.strip())
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"文件不存在: {filepath}")

    if filepath.suffix == ".gz":
        # gzip 文件 → 以文本模式打开，编码 UTF-8
        fh = gzip.open(filepath, mode=mode, encoding="utf-8")
    else:
        # 普通文本文件
        fh = open(filepath, mode=mode, encoding="utf-8")

    try:
        yield fh
    finally:
        fh.close()


def detect_format(filepath: PathLike) -> FileFormat:
    """
    通过读取文件首个非空字符，自动判断文件格式。

    - 以 '>' 开头 → FASTA
    - 以 '@' 开头 → FASTQ

    Parameters
    ----------
    filepath : str or Path
        序列文件路径。

    Returns
    -------
    FileFormat
        检测到的文件格式枚举值。

    Raises
    ------
    ValueError
        如果无法识别文件格式。
    """
    with open_sequence_file(filepath) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue  # 跳过空行
            if stripped.startswith(">"):
                return FileFormat.FASTA
            elif stripped.startswith("@"):
                return FileFormat.FASTQ
            else:
                break

    raise ValueError(
        f"无法识别文件 '{filepath}' 的格式：首个非空行不以 '>' (FASTA) "
        f"或 '@' (FASTQ) 开头"
    )


def infer_format_from_name(filepath: PathLike) -> FileFormat:
    """
    根据文件扩展名推断格式（备用方案，不做内容检测）。

    Parameters
    ----------
    filepath : str or Path

    Returns
    -------
    FileFormat
    """
    filepath = Path(filepath)

    # 去掉 .gz 后缀后再判断
    stem_suffix = filepath.suffix.lower()
    if stem_suffix == ".gz":
        stem_suffix = Path(filepath.stem).suffix.lower()

    if stem_suffix in {".fasta", ".fa", ".fna", ".ffn", ".faa"}:
        return FileFormat.FASTA
    elif stem_suffix in {".fastq", ".fq"}:
        return FileFormat.FASTQ
    else:
        return FileFormat.UNKNOWN