"""
数据模型层：定义序列记录类与统计结果类。

- SeqRecord:  表示一条 FASTA/FASTQ 序列记录（使用 dataclass, frozen=True 保证不可变）
- StatsResult: 表示一条序列的统计结果
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


# ---------------------------------------------------------------------------
# 标准密码子表 (Standard Genetic Code, NCBI Translation Table 1)
# ---------------------------------------------------------------------------
CODON_TABLE: dict[str, str] = {
    # --- 第一第二位 ---
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# DNA 互补碱基映射 (支持 IUPAC 模糊碱基)
COMPLEMENT_MAP: dict[str, str] = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "a": "t", "t": "a", "c": "g", "g": "c",
    "R": "Y", "Y": "R", "S": "S", "W": "W",
    "K": "M", "M": "K", "B": "V", "V": "B",
    "D": "H", "H": "D", "N": "N",
    "r": "y", "y": "r", "s": "s", "w": "w",
    "k": "m", "m": "k", "b": "v", "v": "b",
    "d": "h", "h": "d", "n": "n",
}

# DNA ↔ RNA 转换映射
DNA_TO_RNA: dict[str, str] = {"T": "U", "t": "u"}
RNA_TO_DNA: dict[str, str] = {"U": "T", "u": "t"}


# ---------------------------------------------------------------------------
# 序列记录
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class SeqRecord:
    """
    表示一条 FASTA 或 FASTQ 序列记录。

    Attributes
    ----------
    id : str
        序列标识符（FASTA '>' 后或 FASTQ '@' 后的第一行，不含前缀符号）。
    description : str
        FASTA/FASTQ 头部中除 ID 之外的描述信息。
    seq : str
        核苷酸序列字符串（已去除换行，保持原始大小写）。
    quality : str or None
        FASTQ 质量字符串（Phred+33 编码）。FASTA 记录此字段为 None。
    """

    id: str
    seq: str
    description: str = ""
    quality: Optional[str] = None

    def __post_init__(self) -> None:
        """基本校验。"""
        if self.quality is not None and len(self.quality) != len(self.seq):
            raise ValueError(
                f"SeqRecord '{self.id}': quality length ({len(self.quality)}) "
                f"!= sequence length ({len(self.seq)})"
            )

    @property
    def is_fastq(self) -> bool:
        """是否为 FASTQ 记录。"""
        return self.quality is not None

    def __len__(self) -> int:
        return len(self.seq)

    def to_fasta(self, line_width: int = 80) -> str:
        """将记录格式化为 FASTA 字符串。"""
        lines = [f">{self.id}"]
        seq = self.seq
        for i in range(0, len(seq), line_width):
            lines.append(seq[i : i + line_width])
        return "\n".join(lines)

    def to_fastq(self) -> str:
        """将记录格式化为 FASTQ 字符串（4 行）。"""
        if self.quality is None:
            raise ValueError(f"SeqRecord '{self.id}' 无质量信息，无法输出 FASTQ 格式")
        return f"@{self.id}\n{self.seq}\n+\n{self.quality}"


# ---------------------------------------------------------------------------
# 统计结果
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class StatsResult:
    """
    一条序列的统计结果。

    Attributes
    ----------
    id : str
        序列标识符。
    length : int
        序列总长度。
    gc_count : int
        G + C 碱基总数。
    gc_content : float
        GC 含量（0.0 ~ 1.0）。
    n_count : int
        N（模糊碱基）的数量。
    n_ratio : float
        N 比例（0.0 ~ 1.0）。
    base_counts : dict[str, int]
        各碱基的频数（大写键）。
    """

    id: str
    length: int
    gc_count: int
    gc_content: float
    n_count: int
    n_ratio: float
    base_counts: dict[str, int] = field(default_factory=dict)

    def to_dict(self) -> dict:
        """转为普通字典，方便 CLI 输出 / DataFrame 构建。"""
        return {
            "id": self.id,
            "length": self.length,
            "gc_count": self.gc_count,
            "gc_content": round(self.gc_content, 6),
            "n_count": self.n_count,
            "n_ratio": round(self.n_ratio, 6),
            **{f"count_{base}": cnt for base, cnt in sorted(self.base_counts.items())},
        }

    def summary_line(self) -> str:
        """单行摘要，用于 CLI 快速输出。"""
        return (
            f"{self.id}\tlen={self.length}\t"
            f"GC={self.gc_content:.4f}\t"
            f"N={self.n_ratio:.4f}"
        )