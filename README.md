# bioseqkit
A lightweight, modular Python toolkit for FASTA/FASTQ manipulation built from scratch for bioinformatics learning and lightweight sequence analysis.

## Project Overview
High-throughput sequencing generates massive nucleic acid datasets stored in FASTA and FASTQ formats. Popular bioinformatics libraries like Biopython and SeqKit provide mature parsing APIs but hide underlying low-level logic, making it hard for beginners to understand format specifications, state-machine parsing, and stream-based memory-efficient processing.

This project implements `bioseqkit` entirely from scratch without relying on third-party sequence parsers. It leverages Python generators to achieve constant-memory streaming for large uncompressed and gzipped sequencing files, and encapsulates core bioinformatics algorithms including reverse complement, DNA-RNA transcription, six-frame translation, sequence statistics and k-mer counting. The toolkit ships with a CLI interface, pytest test suites, and Jupyter Notebook visualization pipelines for real genomic data.

This repository serves dual purposes:
1. **Educational training**: Demystify FASTA/FASTQ parsing logic and streaming big-data programming.
2. **Lightweight bioinformatics utility**: Standalone tool for small-to-medium nucleic acid preprocessing and statistical analysis.

## Key Features
### 1. Handwritten FASTA / FASTQ Parser
- Pure self-implemented state-machine parser, no dependency on `Bio.SeqIO`
- Support plain `.fasta`/`.fastq` and gzip-compressed `.fasta.gz`/`.fastq.gz`
- Stream records via Python generators, constant memory footprint for multi-GB files
- Robust fault tolerance: handle multi-line sequences, blank lines, truncated files, ambiguous base `N`, malformed FASTQ quality lines

### 2. Core Sequence Operations
- Reverse complement for DNA strands
- Bidirectional DNA ↔ RNA conversion
- Standard codon table six-frame translation (3 forward frames + 3 reverse-complement frames)
- Mark stop codons with `*`; auto-truncate incomplete triplets at sequence ends

### 3. Sequence Statistics & k-mer Analysis
- Per-record metrics: sequence length, GC content, N ratio, base composition
- Sliding-window k-mer frequency counting
- Output top-N most frequent k-mers
- Lightweight in-memory counter for small genomes; reserved chunked counting interface for ultra-large datasets

### 4. Command-Line Interface (CLI)
Built with `argparse`, subcommands include:
- `stats`: Calculate base statistics and output tabular results
- `revcomp`: Generate reverse-complemented sequences
- `translate`: Perform six-frame protein translation
- Support standard input / standard output pipeline, uniform FASTA output for sequence processing commands

### 5. Testing & Visualization
- Complete pytest test cases covering edge cases: empty files, incomplete records, illegal bases, line-break missing at EOF
- Jupyter Notebook workflow with human chrM mitochondrial genome
- Visualization of base distribution pie charts, GC-content histograms and k-mer frequency plots via Matplotlib & Seaborn

## System Architecture
The project adopts a bottom-up layered design with low coupling and high cohesion:
1. **IO Layer**: Unified handler for plain/gzipped files; state-machine line-by-line parsing; generator-based record yielding
2. **Data Model Layer**: Lightweight `dataclass` Record structure storing ID, sequence, and FASTQ quality strings
3. **Algorithm Layer**: Stateless pure functions for conversion, statistics and k-mer counting
4. **User Interface Layer**: Native Python API + multi-subcommand CLI
5. Test & Visualization Layer: Automated unit tests + real-data exploratory analysis notebooks

## Installation
### 1. Clone Repository
```bash
git clone https://github.com/yourname/bioseqkit.git
cd bioseqkit
```
### 2. Install Package
```bash
pip install bioseqkit
```

## Quick Start
### Use as Python API
```python
from bioseqkit.parser import parse_fastq
from bioseqkit.ops import revcomp, calculate_gc

# Stream FASTQ records without loading all data into memory
for record in parse_fastq("sample.fastq.gz"):
    gc = calculate_gc(record.seq)
    rc_seq = revcomp(record.seq)
    print(f"ID: {record.id}, GC: {gc:.2f}, RevComp: {rc_seq[:20]}...")
```

### Run CLI in Terminal
```bash
# Output sequence statistics
bioseqkit stats input.fasta.gz

# Generate reverse complement sequences
bioseqkit revcomp input.fastq -o output_rc.fasta

# Six-frame translation
bioseqkit translate genome.fasta --frame all
```

### Run Unit Tests
```bash
pytest tests/ -v
```

### Visualization Notebook
Launch Jupyter and open `notebooks/genome_analysis.ipynb` to analyze real genomic data and generate statistical plots.

## Project Structure
```
bioseqkit/
├── src/bioseqkit/        # Core source code
│   ├── parser.py         # FASTA/FASTQ streaming parser
│   ├── record.py         # Sequence Record dataclass model
│   ├── ops.py            # Reverse complement & transcription
│   ├── translate.py      # Six-frame translation logic
│   ├── stats.py          # GC, base count, k-mer functions
│   └── cli.py            # Command-line entry
├── tests/                # Pytest edge-case test suites
├── notebooks/            # Jupyter visualization scripts
├── docs/                 # API documentation
└── pyproject.toml        # Package configuration
```

## License
MIT License
