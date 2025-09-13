# Setup Instructions

## Prerequisites Check

Ensure you have:
- [x] Python 3.12 installed (verified in `.python-version`)
- [x] uv installed and working
- [x] Project dependencies installed (`uv sync`)
- [x] VSCode configured with `.venv` interpreter

## Step 1: Download Practice Data

We'll start with E. coli data as it's small and manageable.

```bash
cd data/raw/

# Download E. coli nanopore reads (~500MB)
wget https://nanopore.s3.climb.ac.uk/MAP006-1_2D_pass.fasta.gz

# Alternative if wget not available:
curl -O https://nanopore.s3.climb.ac.uk/MAP006-1_2D_pass.fasta.gz

# Download reference genome for E. coli
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz
mv GCF_000005845.2_ASM584v2_genomic.fna ../references/ecoli_reference.fasta

cd ../..
```

### Optional: Larger Datasets

```bash
# Zymo Mock Community (~2GB) - 10 known organisms
cd data/raw/
wget https://nanopore.s3.climb.ac.uk/Zymo-GridION-EVEN-BB-SN.fq.gz

# SARS-CoV-2 (small, ~100MB)
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR117/019/SRR11777319/SRR11777319.fastq.gz
```

## Step 2: Create Analysis Scripts

### Script 1: Read Statistics (`scripts/read_stats.py`)

```python
#!/usr/bin/env python3
"""
Calculate basic statistics for nanopore reads
"""

import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from pathlib import Path

def calculate_n50(lengths):
    """Calculate N50 metric for read lengths"""
    sorted_lengths = sorted(lengths, reverse=True)
    total = sum(sorted_lengths)
    target = total / 2
    cumsum = 0
    for length in sorted_lengths:
        cumsum += length
        if cumsum >= target:
            return length
    return 0

def analyze_nanopore_reads(fastq_file, output_dir):
    """
    Analyze nanopore FASTQ/FASTA file and generate statistics
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Determine file format
    file_format = "fastq" if "fastq" in str(fastq_file) or "fq" in str(fastq_file) else "fasta"
    
    # Open file (handle gzip if needed)
    if str(fastq_file).endswith('.gz'):
        handle = gzip.open(fastq_file, "rt")
    else:
        handle = open(fastq_file, "r")
    
    # Parse reads
    print(f"Parsing {file_format} file...")
    reads = []
    qualities = []
    
    for record in SeqIO.parse(handle, file_format):
        reads.append(len(record.seq))
        if file_format == "fastq":
            qualities.append(np.mean(record.letter_annotations["phred_quality"]))
    
    handle.close()
    
    # Calculate statistics
    stats = {
        'Total reads': len(reads),
        'Total bases': sum(reads),
        'Mean length': np.mean(reads),
        'Median length': np.median(reads),
        'Max length': max(reads),
        'Min length': min(reads),
        'N50': calculate_n50(reads),
        'Reads > 10kb': sum(1 for r in reads if r > 10000),
        'Reads > 50kb': sum(1 for r in reads if r > 50000),
    }
    
    if qualities:
        stats['Mean quality'] = np.mean(qualities)
        stats['Median quality'] = np.median(qualities)
    
    # Print statistics
    print("\n=== Nanopore Read Statistics ===")
    for key, value in stats.items():
        if isinstance(value, float):
            print(f"{key}: {value:,.2f}")
        else:
            print(f"{key}: {value:,}")
    
    # Create visualizations
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Read length distribution
    axes[0, 0].hist(reads, bins=50, edgecolor='black')
    axes[0, 0].set_xlabel('Read Length (bp)')
    axes[0, 0].set_ylabel('Count')
    axes[0, 0].set_title('Read Length Distribution')
    axes[0, 0].axvline(stats['N50'], color='red', linestyle='--', label=f'N50: {stats["N50"]:,}')
    axes[0, 0].legend()
    
    # Log-scale read length
    axes[0, 1].hist(reads, bins=50, edgecolor='black')
    axes[0, 1].set_xlabel('Read Length (bp)')
    axes[0, 1].set_ylabel('Count (log scale)')
    axes[0, 1].set_title('Read Length Distribution (Log Scale)')
    axes[0, 1].set_yscale('log')
    
    # Cumulative length plot
    sorted_reads = sorted(reads, reverse=True)
    cumsum = np.cumsum(sorted_reads)
    axes[1, 0].plot(range(len(cumsum)), cumsum)
    axes[1, 0].set_xlabel('Read Rank')
    axes[1, 0].set_ylabel('Cumulative Bases')
    axes[1, 0].set_title('Cumulative Read Length')
    axes[1, 0].axhline(sum(reads)/2, color='red', linestyle='--', label='50% of bases')
    axes[1, 0].legend()
    
    # Quality distribution (if available)
    if qualities:
        axes[1, 1].hist(qualities, bins=30, edgecolor='black')
        axes[1, 1].set_xlabel('Mean Read Quality (Phred)')
        axes[1, 1].set_ylabel('Count')
        axes[1, 1].set_title('Read Quality Distribution')
    else:
        axes[1, 1].text(0.5, 0.5, 'No quality data available\n(FASTA input)', 
                       ha='center', va='center', transform=axes[1, 1].transAxes)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'read_statistics.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Save statistics to file
    with open(output_dir / 'statistics.txt', 'w') as f:
        for key, value in stats.items():
            if isinstance(value, float):
                f.write(f"{key}: {value:,.2f}\n")
            else:
                f.write(f"{key}: {value:,}\n")
    
    return stats

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python read_stats.py <input_fastq> <output_dir>")
        sys.exit(1)
    
    analyze_nanopore_reads(sys.argv[1], sys.argv[2])
```

### Script 2: Simulate Sequencing (`scripts/simulate_sequencing.py`)

```python
#!/usr/bin/env python3
"""
Simulate real-time nanopore sequencing run
"""

import gzip
import time
import random
from datetime import datetime, timedelta
from Bio import SeqIO
import sys

def simulate_nanopore_run(input_file, duration_hours=2, speed_multiplier=1000):
    """
    Simulate reads coming off a nanopore sequencer in real-time
    
    Args:
        input_file: Path to FASTQ/FASTA file
        duration_hours: Simulated run duration
        speed_multiplier: Speed up factor (1000 = 1000x faster than real)
    """
    
    # Determine file format
    file_format = "fastq" if "fastq" in str(input_file) or "fq" in str(input_file) else "fasta"
    
    # Open and count reads
    print("Loading reads...")
    if str(input_file).endswith('.gz'):
        handle = gzip.open(input_file, "rt")
    else:
        handle = open(input_file, "r")
    
    reads = list(SeqIO.parse(handle, file_format))
    handle.close()
    
    total_reads = len(reads)
    reads_per_second = total_reads / (duration_hours * 3600)
    
    print(f"=== Nanopore Sequencing Simulation ===")
    print(f"Total reads to sequence: {total_reads:,}")
    print(f"Simulated duration: {duration_hours} hours")
    print(f"Speed multiplier: {speed_multiplier}x")
    print(f"Reads per second (simulated): {reads_per_second:.2f}")
    print(f"Starting sequencing run...\n")
    
    start_time = datetime.now()
    simulated_start = datetime.now()
    
    bases_sequenced = 0
    reads_completed = 0
    pores_blocked = 0
    
    try:
        for i, read in enumerate(reads):
            # Simulate read production
            bases_sequenced += len(read.seq)
            reads_completed += 1
            
            # Random pore blocking event (characteristic of nanopore)
            if random.random() < 0.001:
                pores_blocked += 1
                print(f"⚠️  Pore blocked! Total blocks: {pores_blocked}")
            
            # Calculate time elapsed
            current_time = datetime.now()
            elapsed = current_time - start_time
            simulated_elapsed = timedelta(seconds=elapsed.total_seconds() * speed_multiplier)
            
            # Print progress every 100 reads
            if i % 100 == 0:
                throughput = bases_sequenced / (elapsed.total_seconds() + 0.001)
                print(f"[{simulated_elapsed}] Reads: {reads_completed:,} | "
                      f"Bases: {bases_sequenced:,} | "
                      f"Throughput: {throughput * speed_multiplier:.0f} bp/s | "
                      f"Active pores: {48 - pores_blocked}")
            
            # Simulate real-time delay
            delay = (1 / reads_per_second) / speed_multiplier
            time.sleep(delay)
            
            # Random early termination (user stops run)
            if random.random() < 0.0001:
                print("\n🛑 User terminated sequencing run early!")
                break
                
    except KeyboardInterrupt:
        print("\n🛑 Sequencing run interrupted by user!")
    
    # Final statistics
    total_time = datetime.now() - start_time
    simulated_total = timedelta(seconds=total_time.total_seconds() * speed_multiplier)
    
    print(f"\n=== Sequencing Run Complete ===")
    print(f"Total reads sequenced: {reads_completed:,}/{total_reads:,}")
    print(f"Total bases: {bases_sequenced:,}")
    print(f"Simulated run time: {simulated_total}")
    print(f"Actual time: {total_time}")
    print(f"Pores blocked during run: {pores_blocked}")
    print(f"Average read length: {bases_sequenced/reads_completed:.0f} bp")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python simulate_sequencing.py <input_file> [duration_hours] [speed_multiplier]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    duration = float(sys.argv[2]) if len(sys.argv) > 2 else 2
    speed = float(sys.argv[3]) if len(sys.argv) > 3 else 1000
    
    simulate_nanopore_run(input_file, duration, speed)
```

### Script 3: Homopolymer Analysis (`scripts/homopolymer_analysis.py`)

```python
#!/usr/bin/env python3
"""
Analyze homopolymer regions in nanopore data
Nanopore sequencers struggle with runs of the same base (AAAAA, TTTTT, etc.)
"""

import gzip
import re
from collections import defaultdict
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

def find_homopolymers(sequence, min_length=4):
    """Find all homopolymer runs in a sequence"""
    homopolymers = []
    pattern = r'(A{' + str(min_length) + r',}|T{' + str(min_length) + r',}|G{' + str(min_length) + r',}|C{' + str(min_length) + r',})'
    
    for match in re.finditer(pattern, str(sequence)):
        homopolymers.append({
            'base': match.group()[0],
            'length': len(match.group()),
            'position': match.start()
        })
    
    return homopolymers

def analyze_homopolymers(input_file, output_prefix="homopolymer"):
    """
    Analyze homopolymer distributions in nanopore reads
    This is important because homopolymer length estimation is a known weakness
    """
    
    file_format = "fastq" if "fastq" in str(input_file) or "fq" in str(input_file) else "fasta"
    
    if str(input_file).endswith('.gz'):
        handle = gzip.open(input_file, "rt")
    else:
        handle = open(input_file, "r")
    
    print("Analyzing homopolymers in reads...")
    
    homopolymer_counts = defaultdict(lambda: defaultdict(int))
    quality_by_homopolymer = defaultdict(list)
    total_reads = 0
    
    for record in SeqIO.parse(handle, file_format):
        total_reads += 1
        homopolymers = find_homopolymers(record.seq)
        
        for hp in homopolymers:
            homopolymer_counts[hp['base']][hp['length']] += 1
            
            # If FASTQ, analyze quality in homopolymer regions
            if file_format == "fastq":
                start = hp['position']
                end = start + hp['length']
                hp_qualities = record.letter_annotations["phred_quality"][start:end]
                quality_by_homopolymer[hp['base']].extend(hp_qualities)
        
        if total_reads % 1000 == 0:
            print(f"Processed {total_reads} reads...")
    
    handle.close()
    
    # Create visualization
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    bases = ['A', 'T', 'G', 'C']
    colors = ['#2ecc71', '#e74c3c', '#f39c12', '#3498db']
    
    for idx, (base, color) in enumerate(zip(bases, colors)):
        ax = axes[idx // 2, idx % 2]
        
        if base in homopolymer_counts:
            lengths = sorted(homopolymer_counts[base].keys())
            counts = [homopolymer_counts[base][l] for l in lengths]
            
            # Only show up to length 20 for clarity
            max_show = min(20, max(lengths) if lengths else 0)
            lengths_show = [l for l in lengths if l <= max_show]
            counts_show = [homopolymer_counts[base][l] for l in lengths_show]
            
            ax.bar(lengths_show, counts_show, color=color, alpha=0.7, edgecolor='black')
            ax.set_xlabel('Homopolymer Length')
            ax.set_ylabel('Count')
            ax.set_title(f'{base} Homopolymers')
            ax.set_yscale('log')
            
            # Add text with statistics
            total_hp = sum(counts)
            avg_length = sum(l * c for l, c in zip(lengths, counts)) / total_hp if total_hp > 0 else 0
            ax.text(0.7, 0.9, f'Total: {total_hp:,}\nAvg length: {avg_length:.1f}', 
                   transform=ax.transAxes, bbox=dict(boxstyle='round', facecolor='wheat'))
    
    plt.suptitle(f'Homopolymer Distribution in Nanopore Reads (n={total_reads:,})', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_distribution.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # If we have quality data, show quality degradation in homopolymers
    if quality_by_homopolymer:
        fig, ax = plt.subplots(figsize=(8, 6))
        
        bp_data = []
        for base in bases:
            if base in quality_by_homopolymer:
                bp_data.append(quality_by_homopolymer[base])
        
        if bp_data:
            bp = ax.boxplot(bp_data, labels=[b for b in bases if b in quality_by_homopolymer])
            ax.set_xlabel('Base')
            ax.set_ylabel('Quality Score (Phred)')
            ax.set_title('Quality Scores in Homopolymer Regions')
            ax.axhline(y=10, color='r', linestyle='--', alpha=0.5, label='Q10 threshold')
            ax.legend()
            plt.savefig(f'{output_prefix}_quality.png', dpi=300, bbox_inches='tight')
            plt.show()
    
    print(f"\nAnalysis complete! Processed {total_reads:,} reads")
    print("\nHomopolymer summary:")
    for base in bases:
        if base in homopolymer_counts:
            total = sum(homopolymer_counts[base].values())
            max_length = max(homopolymer_counts[base].keys())
            print(f"  {base}: {total:,} homopolymers, max length: {max_length}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python homopolymer_analysis.py <input_file>")
        sys.exit(1)
    
    analyze_homopolymers(sys.argv[1])
```

## Step 3: Run Initial Analysis

Test the scripts with the downloaded data:

```bash
# 1. Basic read statistics
uv run python scripts/read_stats.py data/raw/MAP006-1_2D_pass.fasta.gz results/

# 2. Simulate sequencing (runs for 30 seconds at 10000x speed)
uv run python scripts/simulate_sequencing.py data/raw/MAP006-1_2D_pass.fasta.gz 0.5 10000

# 3. Analyze homopolymers
uv run python scripts/homopolymer_analysis.py data/raw/MAP006-1_2D_pass.fasta.gz
```

## Step 4: Create Jupyter Notebook for Exploration

Create `notebooks/explore_nanopore.ipynb`:

```python
# Cell 1: Setup
import sys
sys.path.append('../scripts')
from pathlib import Path
import pandas as pd
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
%matplotlib inline

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (10, 6)

# Cell 2: Load data
data_file = Path("../data/raw/MAP006-1_2D_pass.fasta.gz")
print(f"Loading {data_file.name}...")

import gzip
with gzip.open(data_file, "rt") as handle:
    records = list(SeqIO.parse(handle, "fasta"))
    
print(f"Loaded {len(records):,} reads")
print(f"First read ID: {records[0].id}")
print(f"First read length: {len(records[0].seq):,} bp")

# Cell 3: Read length distribution
read_lengths = [len(r.seq) for r in records]

plt.figure(figsize=(12, 4))
plt.subplot(1, 2, 1)
plt.hist(read_lengths, bins=50, edgecolor='black', alpha=0.7)
plt.xlabel('Read Length (bp)')
plt.ylabel('Count')
plt.title('Read Length Distribution')

plt.subplot(1, 2, 2)
plt.hist(read_lengths, bins=50, edgecolor='black', alpha=0.7, cumulative=True)
plt.xlabel('Read Length (bp)')
plt.ylabel('Cumulative Count')
plt.title('Cumulative Read Length Distribution')
plt.tight_layout()
plt.show()

# Cell 4: Calculate N50
def calculate_n50(lengths):
    sorted_lengths = sorted(lengths, reverse=True)
    total = sum(sorted_lengths)
    cumsum = 0
    for length in sorted_lengths:
        cumsum += length
        if cumsum >= total / 2:
            return length
    return 0

n50 = calculate_n50(read_lengths)
print(f"N50: {n50:,} bp")
print(f"Mean read length: {np.mean(read_lengths):,.0f} bp")
print(f"Median read length: {np.median(read_lengths):,.0f} bp")
print(f"Longest read: {max(read_lengths):,} bp")

# Cell 5: GC content analysis
gc_contents = []
for record in records[:1000]:  # First 1000 reads for speed
    seq = str(record.seq).upper()
    gc_count = seq.count('G') + seq.count('C')
    gc_content = (gc_count / len(seq)) * 100
    gc_contents.append(gc_content)

plt.figure(figsize=(8, 5))
plt.hist(gc_contents, bins=30, edgecolor='black', alpha=0.7, color='steelblue')
plt.xlabel('GC Content (%)')
plt.ylabel('Count')
plt.title('GC Content Distribution (First 1000 Reads)')
plt.axvline(np.mean(gc_contents), color='red', linestyle='--', 
            label=f'Mean: {np.mean(gc_contents):.1f}%')
plt.legend()
plt.show()
```

## Step 5: Start Jupyter Lab

```bash
# Start Jupyter
uv run jupyter lab

# Navigate to notebooks/explore_nanopore.ipynb
```

## Next Steps

### Additional Analyses to Try

1. **Alignment to Reference**
   ```bash
   # Install minimap2 (if not already installed)
   brew install minimap2  # macOS
   # or
   sudo apt-get install minimap2  # Linux
   
   # Align reads
   minimap2 -ax map-ont data/references/ecoli_reference.fasta data/raw/MAP006-1_2D_pass.fasta.gz > results/alignment.sam
   ```

2. **Error Pattern Analysis**
   - Compare aligned reads to reference
   - Identify systematic error patterns
   - Build ML model for error correction

3. **Real-time Analysis Simulation**
   - Process reads in streaming fashion
   - Update statistics as "sequencing progresses"
   - Simulate decision making during a run

### ML Project Ideas

1. **Homopolymer Length Prediction**
   - Train model to correct homopolymer lengths
   - Use context around homopolymers as features
   
2. **Read Quality Prediction**
   - Predict which reads will align well
   - Features: length, GC content, homopolymer count
   
3. **Species Classification**
   - Use Zymo mock community data
   - Build classifier to identify organism from read

### Troubleshooting

**If downloads fail:**
- Try alternative mirrors or use curl instead of wget
- Check available space (need ~2-3GB free)

**If scripts error:**
- Ensure all dependencies installed: `uv sync`
- Check Python version: `uv run python --version` (should be 3.12.x)

**Memory issues with large files:**
- Process reads in chunks rather than loading all at once
- Use generator patterns in Python

## Resources

- [Oxford Nanopore Tutorials](https://community.nanoporetech.com/)
- [BioPython Documentation](https://biopython.org/)
- [Nanopore Analysis Best Practices](https://timkahlke.github.io/LongRead_tutorials/)

## Support

For issues or questions about:
- **Python/environment setup**: Check uv documentation
- **Bioinformatics concepts**: See BioPython tutorials
- **Nanopore specifics**: Visit ONT Community forums