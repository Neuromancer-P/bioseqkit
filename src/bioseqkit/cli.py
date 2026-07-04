"""
src/genbankx/cli.py
命令行接口模块
"""
import argparse
import sys
from typing import List, Dict, Optional
from pathlib import Path

from bioseqkit.parser import parse_fasta, parse_fastq
from bioseqkit.stats import calculate_sequence_stats
from bioseqkit.operations import reverse_complement, six_frame_translation
from bioseqkit.utils import FileFormat, open_sequence_file, detect_format
from bioseqkit.models import SeqRecord


def cmd_stats(args):
    """执行 stats 子命令"""
    file_path = args.input
    output_format = args.format
    
    try:
        # 检测文件格式
        file_format = detect_format(file_path)
        
        # 解析并统计
        results = []
        with open_sequence_file(file_path) as file_handle:
            if file_format == FileFormat.FASTA:
                parser = parse_fasta(file_handle)
            elif file_format == FileFormat.FASTQ:
                parser = parse_fastq(file_handle)
            else:
                print(f"错误：不支持的文件格式 '{file_format}'", file=sys.stderr)
                sys.exit(1)
            
            for seq_record in parser:
                stats = calculate_sequence_stats(seq_record)
                results.append({
                    'id': seq_record.id,
                    'description': seq_record.description,
                    **stats
                })
        
        # 输出结果
        if output_format == 'table':
            print_table(results)
        elif output_format == 'csv':
            print_csv(results)
            
    except FileNotFoundError:
        print(f"错误：文件 '{file_path}' 不存在", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"错误：{e}", file=sys.stderr)
        sys.exit(1)


def cmd_revcomp(args):
    """执行 revcomp 子命令"""
    file_path = args.input
    output_file = args.output
    
    try:
        file_format = detect_format(file_path)
        
        with open_sequence_file(file_path) as file_handle:
            if file_format == FileFormat.FASTA:
                parser = parse_fasta(file_handle)
            elif file_format == FileFormat.FASTQ:
                parser = parse_fastq(file_handle)
            else:
                print(f"错误：不支持的文件格式 '{file_format}'", file=sys.stderr)
                sys.exit(1)
            
            # 输出到文件或终端
            if output_file:
                with open(output_file, 'w') as out_handle:
                    for seq_record in parser:
                        rc_seq = reverse_complement(seq_record.seq)
                        out_handle.write(f">{seq_record.id}_revcomp {seq_record.description}\n")
                        out_handle.write(f"{rc_seq}\n")
                print(f"✅ 反向互补序列已写入：{output_file}")
            else:
                for seq_record in parser:
                    rc_seq = reverse_complement(seq_record.seq)
                    print(f">{seq_record.id}_revcomp {seq_record.description}")
                    print(f"{rc_seq}")
                    
    except FileNotFoundError:
        print(f"错误：文件 '{file_path}' 不存在", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"错误：{e}", file=sys.stderr)
        sys.exit(1)


def cmd_translate(args):
    """执行 translate 子命令"""
    file_path = args.input
    frame = args.frame
    
    try:
        file_format = detect_format(file_path)
        
        with open_sequence_file(file_path) as file_handle:
            if file_format == FileFormat.FASTA:
                parser = parse_fasta(file_handle)
            elif file_format == FileFormat.FASTQ:
                parser = parse_fastq(file_handle)
            else:
                print(f"错误：不支持的文件格式 '{file_format}'", file=sys.stderr)
                sys.exit(1)
            
            for seq_record in parser:
                if frame is None:
                    # 六框翻译
                    proteins = six_frame_translation(seq_record.seq)
                    for i, protein in enumerate(proteins):
                        frame_name = f"Frame {i+1}" if i < 3 else f"RC Frame {i-2}"
                        print(f">{seq_record.id} {frame_name}")
                        print(protein)
                else:
                    # 单框翻译
                    from bioseqkit.operations import translate
                    protein = translate(seq_record.seq, frame)
                    print(f">{seq_record.id} Frame {frame}")
                    print(protein)
                    
    except FileNotFoundError:
        print(f"错误：文件 '{file_path}' 不存在", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"错误：{e}", file=sys.stderr)
        sys.exit(1)


def print_table(results: List[Dict]):
    """以表格形式输出统计结果"""
    if not results:
        print("无数据")
        return
    
    # 定义列
    columns = ['id', 'length', 'gc_content', 'n_ratio', 'A', 'C', 'G', 'T', 'N']
    
    # 打印表头
    header = " | ".join(f"{col:>12}" for col in columns)
    print(header)
    print("-" * len(header))
    
    # 打印数据
    for row in results:
        values = []
        for col in columns:
            value = row.get(col, '')
            if isinstance(value, float):
                values.append(f"{value:>12.2f}")
            else:
                values.append(f"{str(value):>12}")
        print(" | ".join(values))


def print_csv(results: List[Dict]):
    """以 CSV 形式输出统计结果"""
    if not results:
        return
    
    # 定义列
    columns = ['id', 'description', 'length', 'gc_content', 'n_ratio', 'A', 'C', 'G', 'T', 'N']
    
    # 打印表头
    print(",".join(columns))
    
    # 打印数据
    for row in results:
        values = []
        for col in columns:
            value = row.get(col, '')
            # CSV 转义
            if ',' in str(value) or '"' in str(value):
                value = f'"{str(value).replace(chr(34), chr(34)+chr(34))}"'
            values.append(str(value))
        print(",".join(values))


def main():
    """主入口函数"""
    parser = argparse.ArgumentParser(
        prog='bioseqkit',
        description='生物序列处理工具包'
    )
    
    subparsers = parser.add_subparsers(dest='command', help='可用命令')
    
    # stats 子命令
    stats_parser = subparsers.add_parser('stats', help='序列统计分析')
    stats_parser.add_argument('input', help='输入文件路径 (FASTA/FASTQ)')
    stats_parser.add_argument('--format', choices=['table', 'csv'], default='table',
                              help='输出格式 (默认: table)')
    stats_parser.set_defaults(func=cmd_stats)
    
    # revcomp 子命令
    revcomp_parser = subparsers.add_parser('revcomp', help='反向互补序列')
    revcomp_parser.add_argument('input', help='输入文件路径 (FASTA/FASTQ)')
    revcomp_parser.add_argument('--output', '-o', help='输出文件路径 (默认: 终端)')
    revcomp_parser.set_defaults(func=cmd_revcomp)
    
    # translate 子命令
    translate_parser = subparsers.add_parser('translate', help='六框翻译')
    translate_parser.add_argument('input', help='输入文件路径 (FASTA/FASTQ)')
    translate_parser.add_argument('--frame', type=int, choices=[0, 1, 2],
                                  help='阅读框 (0-2，默认: 六框翻译)')
    translate_parser.set_defaults(func=cmd_translate)
    
    # 解析参数
    args = parser.parse_args()
    
    if args.command is None:
        parser.print_help()
        sys.exit(1)
    
    # 执行命令
    args.func(args)


if __name__ == '__main__':
    main()