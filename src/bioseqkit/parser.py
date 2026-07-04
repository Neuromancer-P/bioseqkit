"""
src/genbankx/parser.py
核心序列解析器：支持 FASTA/FASTQ 格式的流式解析
"""
from typing import Generator, Tuple, IO, Optional
from models import SeqRecord
from utils import open_sequence_file, detect_format

def parse_fasta(file_handle: IO[str]) -> Generator[SeqRecord, None, None]:
    """
    解析 FASTA 格式文件，逐条生成 SeqRecord 对象
    
    FASTA 格式示例：
    >seq1 描述信息
    ATCG
    ATCG
    
    Args:
        file_handle: 已打开的文件对象（支持 .fasta / .fa / .fa.gz）
        
    Yields:
        SeqRecord 对象
        
    Raises:
        ValueError: 格式错误时抛出
    """
    header = None
    seq_lines = []  # 累积多行序列
    
    for line in file_handle:
        line = line.strip()
        
        # 跳过空行
        if not line:
            continue
            
        if line.startswith('>'):
            # 遇到新序列：先 yield 上一条
            if header is not None:
                yield SeqRecord(
                    id=header.split()[0],
                    description=' '.join(header.split()[1:]),
                    seq=''.join(seq_lines)
                )
                
            # 重置状态
            header = line[1:].strip()
            seq_lines = []
            
        else:
            # 累积序列行
            if header is None:
                raise ValueError(f"FASTA 格式错误：序列 '{line}' 缺少头部 '>'")
            seq_lines.append(line)
    
    # ⚠️ 关键：yield 最后一条序列（文件结束时的边界处理）
    if header is not None:
        yield SeqRecord(
            id=header.split()[0],
            description=' '.join(header.split()[1:]),
            seq=''.join(seq_lines)
        )


def parse_fastq(file_handle: IO[str]) -> Generator[SeqRecord, None, None]:
    """
    解析 FASTQ 格式文件，逐条生成 SeqRecord 对象
    
    FASTQ 格式示例（每条记录严格4行）：
    @seq1 描述信息
    ATCG
    +
    !!!!
    
    Args:
        file_handle: 已打开的文件对象（支持 .fastq / .fq / .fq.gz）
        
    Yields:
        SeqRecord 对象
        
    Raises:
        ValueError: 格式错误时抛出
    """
    while True:
        # 第1行：header（以 '@' 开头）
        header_line = file_handle.readline().strip()
        
        # 文件结束
        if not header_line:
            break
            
        # 跳过空行（某些 FASTQ 文件记录间有空行）
        if not header_line:
            continue
            
        if not header_line.startswith('@'):
            raise ValueError(
                f"FASTQ 格式错误：期望 '@' 开头，实际得到 '{header_line}'"
            )
        
        # 第2行：序列
        seq_line = file_handle.readline().strip()
        if not seq_line:
            raise ValueError(f"FASTQ 格式错误：header '{header_line}' 后缺少序列")
            
        # 第3行：分隔符（以 '+' 开头，可含可选标识符）
        plus_line = file_handle.readline().strip()
        if not plus_line.startswith('+'):
            raise ValueError(
                f"FASTQ 格式错误：期望 '+' 分隔符，实际得到 '{plus_line}'"
            )
            
        # 第4行：质量值（长度必须等于序列长度）
        qual_line = file_handle.readline().strip()
        if len(qual_line) != len(seq_line):
            raise ValueError(
                f"FASTQ 格式错误：序列长度 {len(seq_line)} ≠ 质量值长度 {len(qual_line)}"
            )
        
        # yield 结果
        yield SeqRecord(
            id=header_line[1:].split()[0],
            description=' '.join(header_line[1:].split()[1:]),
            seq=seq_line,
            qual=qual_line
        )