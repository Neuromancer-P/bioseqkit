from pathlib import Path
from setuptools import setup, find_packages

here = Path(__file__).resolve().parent
long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
    name="bioseqkits",
    version="0.1.0",
    description="Lightweight Python toolkit for FASTA/FASTQ parsing and bioinformatics sequence analysis.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Ouyang Qinglang, Lin Yushu, Liu Yijun",
    author_email="ouyang725@sjtu.edu.cn, benzi0228@sjtu.edu.cn, liu_yijun@sjtu.edu.cn",
    license="MIT",
    python_requires=">=3.9",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[],
    extras_require={
        "test": ["pytest>=7.0"],
        "notebook": ["matplotlib>=3.5", "seaborn>=0.11"],
        "dev": ["pytest>=7.0", "jupyter", "black", "isort"],
    },
    entry_points={
        "console_scripts": [
            "bioseqkits=bioseqkits.cli:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3.10",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
    ],
)
