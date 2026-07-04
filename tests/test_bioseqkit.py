from __future__ import annotations

import sys
from pathlib import Path

project_root = Path(__file__).resolve().parents[1]
module_path = project_root / "src" / "bioseqkit"
sys.path.insert(0, str(module_path))

from utils import open_sequence_file, detect_format, FileFormat
from parser import parse_fasta
from stats import calculate_sequence_stats
from operations import reverse_complement, translate, six_frame_translation


def test_import_modules() -> None:
    assert callable(open_sequence_file)
    assert callable(parse_fasta)
    assert callable(calculate_sequence_stats)
    assert callable(reverse_complement)
    assert callable(translate)
    assert callable(six_frame_translation)


def test_detect_format_fasta(sequence_fasta: Path) -> None:
    assert detect_format(sequence_fasta) == FileFormat.FASTA


def test_parse_fasta_record(sequence_fasta: Path) -> None:
    with open_sequence_file(sequence_fasta) as fh:
        records = list(parse_fasta(fh))

    assert len(records) == 1
    record = records[0]
    assert record.id == "NC_012920.1"
    assert record.seq.startswith("GATCACAGGT")
    assert len(record.seq) > 100


def test_calculate_sequence_stats(sequence_fasta: Path) -> None:
    with open_sequence_file(sequence_fasta) as fh:
        record = next(parse_fasta(fh))

    stats = calculate_sequence_stats(record)
    assert stats["length"] == len(record.seq)
    assert 0 <= stats["gc_content"] <= 100
    assert 0 <= stats["n_ratio"] <= 100
    assert set(stats) >= {"length", "gc_content", "n_ratio", "A", "C", "G", "T", "N"}


def test_reverse_complement_simple() -> None:
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("atgc") == "gcat"


def test_translate_simple() -> None:
    assert translate("ATGGCC", 0) == "MA"
    assert translate("ATGGCC", 1) == "W"
    assert translate("ATGGCC", 2) == "G"


def test_six_frame_translation_simple() -> None:
    result = six_frame_translation("ATGGCC")
    assert len(result) == 6
    assert result[0] == "MA"
    assert result[3] == "GG"
