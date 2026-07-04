"""
src/genbankx/operations.py
序列操作算法模块
"""
from models import SeqRecord
from typing import Dict, List, Tuple, Optional

# 碱基互补映射
COMPLEMENT_MAP = {
    'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
    'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
    'N': 'N', 'n': 'n'
}

# RNA 与 DNA 互转映射
DNA_TO_RNA = {'T': 'U', 't': 'u'}
RNA_TO_DNA = {'U': 'T', 'u': 't'}

# 遗传密码表 (标准遗传密码)
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def reverse_complement(seq: str) -> str:
    """计算序列的反向互补链"""
    return ''.join(COMPLEMENT_MAP.get(base, base) for base in reversed(seq))

def dna_to_rna(seq: str) -> str:
    """将 DNA 序列转换为 RNA 序列"""
    return ''.join(DNA_TO_RNA.get(base, base) for base in seq)

def rna_to_dna(seq: str) -> str:
    """将 RNA 序列转换为 DNA 序列"""
    return ''.join(RNA_TO_DNA.get(base, base) for base in seq)

def translate(seq: str, frame: int = 0) -> str:
    """
    将 DNA 序列翻译为氨基酸序列
    
    Args:
        seq: DNA 序列
        frame: 阅读框 (0, 1, 2)
        
    Returns:
        氨基酸序列
    """
    # 确保序列长度能被3整除
    seq = seq[frame:]
    seq = seq[:len(seq) // 3 * 3]
    
    # 翻译
    protein = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3].upper()
        amino_acid = CODON_TABLE.get(codon, 'X')
        protein.append(amino_acid)
        
    return ''.join(protein)

def six_frame_translation(seq: str) -> List[str]:
    """
    六框翻译：正向3框 + 反向互补3框
    
    Args:
        seq: DNA 序列
        
    Returns:
        6 个阅读框的翻译结果
    """
    rc_seq = reverse_complement(seq)
    return [
        translate(seq, 0),
        translate(seq, 1),
        translate(seq, 2),
        translate(rc_seq, 0),
        translate(rc_seq, 1),
        translate(rc_seq, 2)
    ]

def kmer_frequency(seq: str, k: int = 3) -> Dict[str, int]:
    """
    计算 k-mer 频率
    
    Args:
        seq: 序列
        k: k-mer 长度
        
    Returns:
        k-mer 频率字典
    """
    kmer_counts = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
    return kmer_counts

def top_k_kmers(kmer_freq: Dict[str, int], top_k: int = 10) -> List[Tuple[str, int]]:
    """
    获取 top-k k-mer
    
    Args:
        kmer_freq: k-mer 频率字典
        top_k: 返回前 top_k 个
        
    Returns:
        (kmer, count) 元组列表
    """
    return sorted(kmer_freq.items(), key=lambda x: x[1], reverse=True)[:top_k]