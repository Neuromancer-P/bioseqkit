# bioseqkit

A lightweight modular Python toolkit for FASTA/FASTQ parsing and sequence analysis, implemented from scratch with pure Python.

## Project Overview

`bioseqkit` is a teaching-oriented toolkit that recreates core FASTA/FASTQ parsing and sequence analysis logic without relying on `Bio.SeqIO` or other sequence parsing libraries. It uses Python generators to parse sequence files in a streaming fashion, supports plain and gzip-compressed input, and exposes common bioinformatics operations through both a Python API and a command-line interface.

This repository currently includes:
- `src/bioseqkit/` package with parser, sequence model, operations, stats, and CLI modules
- pytest-based unit tests under `tests/`
- Jupyter notebooks for interactive sequence analysis and validation

## Key Features

### 1. FASTA / FASTQ Parsing
- Pure Python parser for FASTA and FASTQ files
- Supports both plain and gzip-compressed files
- Generator-based streaming parsing for low memory usage
- Handles multi-line sequences, blank lines, and common FASTQ formatting edge cases

### 2. Sequence Operations
- `reverse_complement` for DNA reverse complement computation
- `dna_to_rna` and `rna_to_dna` conversion
- `translate` for frame-specific translation
- `six_frame_translation` for all six reading frames
- `kmer_frequency` and `top_k_kmers` for k-mer analysis

### 3. Sequence Statistics
- Per-record sequence length
- GC content calculation
- N base ratio calculation
- Base composition percentages for A/C/G/T/N
- Convenience functions like `calculate_sequence_stats`

### 4. Command-Line Interface (CLI)
Supported subcommands:
- `stats` — sequence statistics table output
- `revcomp` — reverse complement sequence generation
- `translate` — protein translation of input sequences

### 5. Testing and Notebooks
- `pytest` test suite for parser and algorithm edge cases
- Notebook examples for FASTA analysis and validation
- Current notebook file: `analysis_sequence_fasta.ipynb`

## Installation

### Clone repository
```bash
git clone https://github.com/Neuromancer-P/bioseqkit.git
cd bioseqkit
```

### Install package
```bash
pip install .
```

### Install development dependencies
```bash
pip install -e .[dev]
```

## Quick Start

### Python API example
```python
from bioseqkit.parser import parse_fasta
from bioseqkit.operations import reverse_complement, dna_to_rna
from bioseqkit.stats import calculate_sequence_stats

with open('tests/data/sequence.fasta') as fh:
    for record in parse_fasta(fh):
        stats = calculate_sequence_stats(record)
        rc = reverse_complement(record.seq)
        rna = dna_to_rna(record.seq)
        print(record.id, stats['length'], stats['gc_content'], rc[:20], rna[:20])
```

### Run CLI commands
```bash
bioseqkit stats tests/data/sequence.fasta
bioseqkit revcomp tests/data/sequence.fasta -o output_revcomp.fasta
bioseqkit translate tests/data/sequence.fasta --frame 0
```

### Run tests
```bash
pytest tests/ -v
```

## Project Structure

```
src/bioseqkit/
├── cli.py
├── models.py
├── operations.py
├── parser.py
├── stats.py
├── utils.py
├── __init__.py
tests/
├── __init__.py
├── conftest.py
├── test_bioseqkit.py
├── data/
│   └── sequence.fasta
analysis_sequence_fasta.ipynb
bioseqkit_test.ipynb
pyproject.toml
README.md
setup.py
```

## Packaging

This project is configured with `pyproject.toml` and uses `setuptools` for packaging. The package exposes a console script named `bioseqkit` that points to `bioseqkit.cli:main`.

## License

MIT License
