"""
src/genbankx/stats.py
序列统计算法模块
"""
from bioseqkits.models import SeqRecord
from bioseqkits.utils import FileFormat, open_sequence_file, detect_format
from typing import Dict, List

def calculate_length(seq: str) -> int:
    """计算序列长度"""
    return len(seq)

def calculate_gc_content(seq: str) -> float:
    """计算 GC 含量 (0-100%)"""
    gc_count = sum(1 for base in seq.upper() if base in 'GC')
    return (gc_count / len(seq)) * 100 if len(seq) > 0 else 0.0

def calculate_n_ratio(seq: str) -> float:
    """计算 N 比例 (0-100%)"""
    n_count = sum(1 for base in seq.upper() if base == 'N')
    return (n_count / len(seq)) * 100 if len(seq) > 0 else 0.0

def calculate_base_composition(seq: str) -> Dict[str, float]:
    """计算碱基组成比例"""
    seq_upper = seq.upper()
    total = len(seq_upper)
    
    if total == 0:
        return {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0, 'N': 0.0}
    
    return {
        'A': (seq_upper.count('A') / total) * 100,
        'C': (seq_upper.count('C') / total) * 100,
        'G': (seq_upper.count('G') / total) * 100,
        'T': (seq_upper.count('T') / total) * 100,
        'N': (seq_upper.count('N') / total) * 100
    }

def calculate_sequence_stats(seq_record: SeqRecord) -> Dict[str, float]:
    """计算单条序列的完整统计信息"""
    return {
        'length': calculate_length(seq_record.seq),
        'gc_content': calculate_gc_content(seq_record.seq),
        'n_ratio': calculate_n_ratio(seq_record.seq),
        **calculate_base_composition(seq_record.seq)
    }

def calculate_stats_from_file(file_path: str) -> List[Dict]:
    """
    从文件计算所有序列的统计信息
    
    Args:
        file_path: FASTA/FASTQ 文件路径
        
    Returns:
        每条序列的统计信息列表
    """
    from bioseqkits.parser import parse_fasta, parse_fastq
    
    stats_list = []
    file_format = detect_format(file_path)
    
    with open_sequence_file(file_path) as file_handle:
        if file_format == FileFormat.FASTA:
            for seq_record in parse_fasta(file_handle):
                stats_list.append({
                    'id': seq_record.id,
                    'description': seq_record.description,
                    **calculate_sequence_stats(seq_record)
                })
        elif file_format == FileFormat.FASTQ:
            for seq_record in parse_fastq(file_handle):
                stats_list.append({
                    'id': seq_record.id,
                    'description': seq_record.description,
                    **calculate_sequence_stats(seq_record)
                })
    
    return stats_list