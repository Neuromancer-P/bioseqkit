# tests/conftest.py
"""pytest 全局 fixtures:提供临时测试数据文件。"""

import tempfile
from pathlib import Path

import pytest


@pytest.fixture()
def tmp_dir(tmp_path: Path) -> Path:
    """返回 pytest 提供的临时目录。"""
    return tmp_path


@pytest.fixture()
def simple_fasta(tmp_dir: Path) -> Path:
    """创建一个简单的多行 FASTA 文件。"""
    p = tmp_dir / "simple.fasta"
    p.write_text(
        ">seq1 description here\n"
        "ATCGATCG\n"
        "GGCCTTAA\n"
        ">seq2\n"
        "NNNNATCG\n"
        "\n"  # 空行（边界情况）
        ">seq3\n"
        "ATCG\n",
        encoding="utf-8",
    )
    return p


@pytest.fixture()
def simple_fastq(tmp_dir: Path) -> Path:
    """创建一个简单的 FASTQ 文件。"""
    p = tmp_dir / "simple.fastq"
    p.write_text(
        "@read1\n"
        "ATCGATCG\n"
        "+\n"
        "IIIIIIII\n"
        "@read2\n"
        "GGCCNNAA\n"
        "+\n"
        "IIIIIIII\n",
        encoding="utf-8",
    )
    return p


@pytest.fixture()
def gzipped_fasta(simple_fasta: Path, tmp_dir: Path) -> Path:
    """创建 simple_fasta 的 gzip 压缩版本。"""
    import gzip

    gz_path = tmp_dir / "simple.fasta.gz"
    with open(simple_fasta, "rb") as f_in:
        with gzip.open(gz_path, "wb") as f_out:
            f_out.write(f_in.read())
    return gz_path


@pytest.fixture()
def sequence_fasta() -> Path:
    """返回 tests/data/sequence.fasta 的路径。"""
    return Path(__file__).parent / "data" / "sequence.fasta"