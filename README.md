# Nanopore Data Analysis

A Python project for exploring and analyzing Oxford Nanopore sequencing data without requiring actual sequencing hardware.

## Overview

This project provides tools and workflows for:
- Analyzing nanopore read statistics and quality metrics
- Simulating real-time sequencing runs
- Identifying and analyzing homopolymer regions (a known nanopore challenge)
- Working with real publicly available nanopore datasets

## Quick Start

```bash
# Install dependencies
uv sync --group dev

# Run basic analysis
uv run python scripts/analyze_reads.py data/raw/ecoli_nanopore_reads.fastq.gz results/

# Start Jupyter for interactive analysis
uv run jupyter lab

# Code quality (linting & formatting)
uv run ruff check          # Check for linting issues
uv run ruff check --fix    # Auto-fix issues where possible
uv run ruff format         # Format all code
```

## Project Structure

```
nanopore/
├── data/              # Sequencing data
│   ├── raw/          # Original downloaded files
│   ├── processed/    # Processed FASTQ/FASTA files
│   └── references/   # Reference genomes
├── notebooks/        # Jupyter notebooks for exploration
├── scripts/          # Analysis scripts
├── results/          # Analysis outputs
│   ├── figures/     # Visualizations
│   └── reports/     # Analysis reports
└── README.md        # This file
```

## Requirements

- Python 3.12
- ~10GB free disk space for data
- macOS/Linux (Windows users should use WSL2)

## Setup

See [notes/setup_instructions.md](notes/setup_instructions.md) for detailed setup and data download instructions.

### Code Quality

This project uses **Ruff** for fast linting and formatting, with **pre-commit hooks** for automatic quality checks:

```bash
# Manual linting & formatting
uv run ruff check --fix    # Lint and auto-fix issues
uv run ruff format         # Format all Python code and notebooks

# Pre-commit setup (runs automatically on git commits)
uv run pre-commit install  # Install hooks
uv run pre-commit run --all-files  # Run on all files manually
```

**VS Code Integration**: Install the [Ruff extension](https://marketplace.visualstudio.com/items?itemName=charliermarsh.ruff) for real-time linting and formatting.

## Key Features

- **No Hardware Required**: Work with real nanopore data from public repositories
- **ML-Ready**: Designed for ML engineers new to bioinformatics
- **Real Data**: Use actual sequencing data with authentic error profiles and characteristics
- **Educational**: Learn nanopore-specific challenges like homopolymer errors and long-read analysis

## Technologies

- **uv**: Fast Python package management
- **BioPython**: Sequence analysis and file parsing
- **ONT Tools**: Official Oxford Nanopore Python libraries
- **Jupyter**: Interactive data exploration
- **Ruff**: Lightning-fast Python linting and formatting
- **Pre-commit**: Automated code quality checks

## License

MIT
