# Nanopore Data Analysis

A hands-on learning project for understanding Oxford Nanopore sequencing data analysis - designed for ML engineers new to bioinformatics.

## What is Nanopore Sequencing?

**New to nanopore sequencing?** Here's the key idea: Imagine DNA passing through a tiny hole (nanopore) one nucleotide at a time. Each base (A, T, G, C) creates a unique electrical signal that we can detect and decode into sequence data.

**Why nanopore is exciting:**
- **Long reads**: 1kb to 100kb+ sequences (vs ~150bp for traditional sequencing)
- **Real-time**: Get data as it sequences (no waiting for batch processing)
- **Direct RNA**: Can sequence RNA directly without conversion to DNA

**Why nanopore is challenging:**
- **Higher error rates**: ~5-15% vs ~0.1% for short-read sequencing
- **Homopolymer errors**: Struggles with repetitive bases like "AAAAA"
- **Complex data**: Raw signals require sophisticated processing

## What This Project Does

This is a **hands-on learning tool** for understanding nanopore data analysis:

- **📊 Analyze real data**: Work with authentic E. coli nanopore sequences
- **🧬 Learn error patterns**: See homopolymer errors and quality issues firsthand
- **📈 Understand metrics**: Length distributions, quality scores, N50 statistics
- **🎓 ML-friendly**: Explanations designed for engineers new to bioinformatics

**What we DON'T do:** Raw signal processing (FAST5 → FASTQ). We start with pre-processed FASTQ files to focus on the biological insights rather than signal processing complexity.

## 🚀 Getting Started

### 1. Setup Environment
```bash
# Clone and setup
git clone <your-repo>
cd nanopore

# Install dependencies (including dev tools)
uv sync --group dev

# Setup code quality hooks (optional but recommended)
uv run pre-commit install
```

### 2. Download Data
```bash
# Download real E. coli nanopore data (~2GB)
uv run python scripts/download_data.py

# Verify the download worked
uv run python scripts/verify_data.py
```

### 3. Start Analyzing!

**Option A: Interactive Notebook** (Recommended for learning)
```bash
# Launch Jupyter and open the exploration notebook
uv run jupyter lab notebooks/
```
The notebook walks you through:
- Understanding nanopore data structure
- Quality analysis and visualization
- Identifying homopolymer errors
- Statistical analysis with biological context

**Option B: Command Line Analysis**
```bash
# Run automated analysis pipeline
uv run python scripts/analyze_reads.py data/raw/ecoli_nanopore_reads.fastq.gz results/
```

### 4. Code Quality (Optional)
```bash
# Format and lint your code
uv run ruff check --fix    # Auto-fix linting issues
uv run ruff format         # Format all Python code

# Or just commit - pre-commit hooks will run automatically!
git add . && git commit -m "your changes"
```

## 📁 Project Structure

```
nanopore/
├── 📊 data/              # Sequencing data (created after download)
│   ├── raw/             # Original FASTQ files (~2GB)
│   └── references/      # E. coli reference genome
├── 📓 notebooks/        # Interactive analysis
│   └── nanopore_exploration.ipynb  # Main tutorial notebook
├── 🐍 scripts/          # Analysis tools
│   ├── download_data.py    # Get real nanopore data
│   ├── verify_data.py      # Check data integrity
│   └── analyze_reads.py    # Automated analysis pipeline
├── 📈 results/          # Analysis outputs (created during analysis)
│   ├── figures/        # Plots and visualizations
│   └── reports/        # Statistical summaries
└── ⚙️  config files     # Development setup
    ├── pyproject.toml     # Dependencies & ruff config
    ├── .pre-commit-config.yaml  # Code quality hooks
    └── .vscode/settings.json    # VS Code integration
```

## Requirements

- Python 3.12
- ~10GB free disk space for data
- macOS/Linux (Windows users should use WSL2)

## 📚 What You'll Learn

### Core Nanopore Concepts
- **Read length distributions**: Why nanopore produces variable-length sequences
- **Quality scores**: How to interpret nanopore quality metrics (different from Illumina!)
- **Error patterns**: Why "AAAAA" sequences cause problems (homopolymer errors)
- **N50 statistics**: The key metric for long-read sequencing quality

### Practical Skills
- **FASTQ file handling**: Working with large bioinformatics file formats
- **Statistical analysis**: Length distributions, quality correlations, error rates
- **Data visualization**: Creating publication-quality plots of genomic data
- **Python + Biology**: Using pandas, matplotlib, and BioPython for genomics

### Why Start with FASTQ?
We deliberately start with **processed FASTQ files** instead of raw FAST5 signal data because:
- **Focus on biology**: Learn sequence analysis without signal processing complexity
- **Standard format**: FASTQ is the universal format for downstream analysis
- **Manageable size**: FASTQ files are smaller and easier to work with than raw signals
- **Real-world workflow**: Most nanopore analysis starts with basecalled FASTQ data

## 🛠️ Advanced Setup

### Development Environment
```bash
# VS Code users: Install the Ruff extension for real-time linting
# https://marketplace.visualstudio.com/items?itemName=charliermarsh.ruff

# Full development setup
uv run pre-commit install              # Auto-format on commit
uv run pre-commit run --all-files      # Format everything now
```

### Troubleshooting
- **Large files**: The data download is ~2GB. Ensure you have space and good internet.
- **Windows**: Use WSL2 for best compatibility with bioinformatics tools.
- **Memory**: Large datasets may need 8GB+ RAM for full analysis.

**Having issues?** Check the troubleshooting section above or open an issue.

## ✨ Key Features

- **🧬 Real Data**: Work with authentic E. coli nanopore sequencing data
- **🎓 ML-Friendly**: Designed for ML engineers learning bioinformatics
- **📊 Interactive**: Jupyter notebooks with rich visualizations
- **⚡ Modern Tools**: uv + Ruff for fast development workflow
- **🔍 Educational**: Learn nanopore-specific error patterns and analysis techniques
- **💻 No Hardware**: No sequencing equipment needed - just Python!

## Technologies

- **uv**: Fast Python package management
- **BioPython**: Sequence analysis and file parsing
- **ONT Tools**: Official Oxford Nanopore Python libraries
- **Jupyter**: Interactive data exploration
- **Ruff**: Lightning-fast Python linting and formatting
- **Pre-commit**: Automated code quality checks

## License

MIT
