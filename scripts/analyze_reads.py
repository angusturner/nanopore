#!/usr/bin/env python3
"""
Analyze nanopore sequencing reads - Educational version for ML engineers.

This script provides a comprehensive analysis of nanopore reads, explaining
each concept from both bioinformatics and ML perspectives.
"""

import gzip
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from Bio import SeqIO

warnings.filterwarnings("ignore")

# Set up plotting
plt.style.use("default")
sns.set_palette("husl")


def load_reads(fastq_file):
    """
    Load nanopore reads from FASTQ/FASTA file.

    For ML engineers:
    - Think of each read as a data sample
    - Each nucleotide (A, T, G, C) is like a categorical feature
    - Read length is a continuous feature
    - Quality scores (if FASTQ) are like confidence scores
    - We're dealing with variable-length sequences (like text/NLP)
    """
    print(f"🧬 Loading reads from {fastq_file}...")

    reads_data = []

    # Determine format from filename
    file_format = "fastq" if "fastq" in str(fastq_file) or "fq" in str(fastq_file) else "fasta"
    print(f"📄 Detected format: {file_format.upper()}")

    # Handle compressed files and parse sequences
    open_func = gzip.open if str(fastq_file).endswith(".gz") else open
    open_mode = "rt" if str(fastq_file).endswith(".gz") else "r"

    with open_func(fastq_file, open_mode) as handle:
        for i, record in enumerate(SeqIO.parse(handle, file_format)):
            read_data = {
                "read_id": record.id,
                "sequence": str(record.seq),
                "length": len(record.seq),
            }

            # Add quality scores if FASTQ
            if (
                file_format == "fastq"
                and hasattr(record, "letter_annotations")
                and "phred_quality" in record.letter_annotations
            ):
                read_data["quality_scores"] = record.letter_annotations["phred_quality"]
                read_data["mean_quality"] = sum(record.letter_annotations["phred_quality"]) / len(
                    record.letter_annotations["phred_quality"]
                )

            reads_data.append(read_data)

            # Progress indicator
            if (i + 1) % 5000 == 0:
                print(f"  Loaded {i + 1:,} reads...")

    print(f"✅ Loaded {len(reads_data):,} total reads")
    has_quality = any("quality_scores" in read for read in reads_data)
    if has_quality:
        print("📊 Quality scores available: ✅ Yes! (FASTQ format)")
    else:
        print("📊 Quality scores available: ❌ No (FASTA format)")

    return reads_data


def calculate_basic_stats(reads_data):
    """
    Calculate basic statistics about the sequencing run.

    ML Context:
    - These are like dataset statistics you'd calculate before training
    - Help understand data distribution and quality
    - N50 is like a median but weighted by length (important for genomics)
    """

    lengths = [read["length"] for read in reads_data]

    # Calculate N50 - a key genomics metric
    def calculate_n50(lengths):
        """
        N50: Length such that 50% of total bases are in reads >= this length

        Think of it as: "Half the sequenced DNA comes from reads this long or longer"
        Different from median because it's weighted by read length.
        """
        sorted_lengths = sorted(lengths, reverse=True)
        total_bases = sum(sorted_lengths)
        cumulative_bases = 0

        for length in sorted_lengths:
            cumulative_bases += length
            if cumulative_bases >= total_bases / 2:
                return length
        return 0

    stats = {
        "total_reads": len(reads_data),
        "total_bases": sum(lengths),
        "mean_length": np.mean(lengths),
        "median_length": np.median(lengths),
        "std_length": np.std(lengths),
        "min_length": min(lengths),
        "max_length": max(lengths),
        "n50": calculate_n50(lengths),
    }

    # Add percentiles
    for p in [5, 25, 75, 95]:
        stats[f"p{p}_length"] = np.percentile(lengths, p)

    return stats


def analyze_sequence_composition(reads_data, sample_size=10000):
    """
    Analyze the nucleotide composition of reads.

    ML Context:
    - Like analyzing token frequencies in NLP
    - GC content is important for genomics (affects DNA stability)
    - Can reveal sequencing biases or contamination
    """

    print(
        f"🔬 Analyzing sequence composition (sampling {min(sample_size, len(reads_data)):,} reads)..."
    )

    # Sample reads for efficiency
    sample_reads = reads_data[:sample_size] if len(reads_data) > sample_size else reads_data

    # Nucleotide counts
    nucleotide_counts = {"A": 0, "T": 0, "G": 0, "C": 0, "N": 0}
    gc_contents = []

    for read in sample_reads:
        seq = read["sequence"].upper()

        # Count nucleotides
        for nuc in nucleotide_counts:
            nucleotide_counts[nuc] += seq.count(nuc)

        # Calculate GC content for this read
        gc_count = seq.count("G") + seq.count("C")
        at_count = seq.count("A") + seq.count("T")
        total_valid = gc_count + at_count

        if total_valid > 0:
            gc_content = (gc_count / total_valid) * 100
            gc_contents.append(gc_content)

    # Calculate overall composition
    total_nucleotides = sum(nucleotide_counts.values())
    composition = {nuc: count / total_nucleotides * 100 for nuc, count in nucleotide_counts.items()}

    return {
        "nucleotide_counts": nucleotide_counts,
        "composition_percent": composition,
        "gc_contents": gc_contents,
        "mean_gc_content": np.mean(gc_contents),
        "std_gc_content": np.std(gc_contents),
    }


def find_homopolymers(sequence, min_length=4):
    """
    Find homopolymer runs (AAAA, TTTT, etc.)

    Why this matters:
    - Nanopore sequencers struggle with repetitive sequences
    - Like how RNNs struggle with long repeated patterns
    - Major source of sequencing errors
    """
    import re

    homopolymers = []
    for base in ["A", "T", "G", "C"]:
        pattern = f"{base}{{{min_length},}}"  # e.g., 'A{4,}' matches AAAA+

        for match in re.finditer(pattern, sequence):
            homopolymers.append(
                {"base": base, "length": match.end() - match.start(), "position": match.start()}
            )

    return homopolymers


def analyze_homopolymers(reads_data, sample_size=5000):
    """
    Analyze homopolymer distribution.

    This is nanopore-specific analysis - these sequencers have trouble
    with runs of identical bases (like AAAAA or TTTTT).
    """

    print(f"🔄 Analyzing homopolymers (sampling {min(sample_size, len(reads_data)):,} reads)...")

    sample_reads = reads_data[:sample_size] if len(reads_data) > sample_size else reads_data

    homopolymer_stats = {"A": [], "T": [], "G": [], "C": []}

    for read in sample_reads:
        homopolymers = find_homopolymers(read["sequence"])
        for hp in homopolymers:
            homopolymer_stats[hp["base"]].append(hp["length"])

    # Calculate statistics for each base
    hp_summary = {}
    for base in ["A", "T", "G", "C"]:
        lengths = homopolymer_stats[base]
        if lengths:
            hp_summary[base] = {
                "count": len(lengths),
                "mean_length": np.mean(lengths),
                "max_length": max(lengths),
                "lengths": lengths,
            }
        else:
            hp_summary[base] = {"count": 0, "mean_length": 0, "max_length": 0, "lengths": []}

    return hp_summary


def create_visualizations(stats, composition, homopolymer_data, reads_data, output_dir="results"):
    """
    Create comprehensive visualizations of the sequencing data.
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    # Create figure with subplots
    plt.figure(figsize=(16, 12))

    # 1. Read length distribution
    plt.subplot(2, 3, 1)
    lengths = [read["length"] for read in reads_data]
    plt.hist(lengths, bins=50, edgecolor="black", alpha=0.7)
    plt.axvline(stats["n50"], color="red", linestyle="--", label=f"N50: {stats['n50']:,.0f} bp")
    plt.axvline(
        stats["mean_length"],
        color="orange",
        linestyle="--",
        label=f"Mean: {stats['mean_length']:,.0f} bp",
    )
    plt.xlabel("Read Length (bp)")
    plt.ylabel("Count")
    plt.title("Read Length Distribution")
    plt.legend()
    plt.yscale("log")

    # 2. Cumulative read length
    plt.subplot(2, 3, 2)
    sorted_lengths = sorted(lengths, reverse=True)
    cumulative = np.cumsum(sorted_lengths)
    plt.plot(range(len(cumulative)), cumulative, linewidth=2)
    plt.axhline(stats["total_bases"] / 2, color="red", linestyle="--", label="50% of bases")
    plt.xlabel("Read Rank")
    plt.ylabel("Cumulative Bases")
    plt.title("Cumulative Length Distribution")
    plt.legend()

    # 3. Nucleotide composition
    plt.subplot(2, 3, 3)
    bases = ["A", "T", "G", "C"]
    percentages = [composition["composition_percent"][base] for base in bases]
    colors = ["#ff9999", "#66b3ff", "#99ff99", "#ffcc99"]

    bars = plt.bar(bases, percentages, color=colors, edgecolor="black")
    plt.ylabel("Percentage")
    plt.title("Nucleotide Composition")
    plt.ylim(0, 35)

    # Add percentage labels on bars
    for bar, pct in zip(bars, percentages, strict=False):
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.5,
            f"{pct:.1f}%",
            ha="center",
            va="bottom",
        )

    # 4. GC content distribution
    plt.subplot(2, 3, 4)
    plt.hist(composition["gc_contents"], bins=30, edgecolor="black", alpha=0.7)
    plt.axvline(
        composition["mean_gc_content"],
        color="red",
        linestyle="--",
        label=f"Mean: {composition['mean_gc_content']:.1f}%",
    )
    plt.xlabel("GC Content (%)")
    plt.ylabel("Count")
    plt.title("GC Content Distribution")
    plt.legend()

    # 5. Homopolymer analysis
    plt.subplot(2, 3, 5)
    hp_counts = [homopolymer_data[base]["count"] for base in bases]
    bars = plt.bar(bases, hp_counts, color=colors, edgecolor="black")
    plt.ylabel("Homopolymer Count")
    plt.title("Homopolymer Frequency (≥4 bases)")
    plt.yscale("log")

    # 6. Read length vs GC content (scatter sample)
    plt.subplot(2, 3, 6)
    sample_size = min(5000, len(reads_data))
    sample_lengths = [reads_data[i]["length"] for i in range(sample_size)]
    sample_gc = composition["gc_contents"][: len(sample_lengths)]

    plt.scatter(sample_gc, sample_lengths, alpha=0.5, s=1)
    plt.xlabel("GC Content (%)")
    plt.ylabel("Read Length (bp)")
    plt.title(f"Length vs GC Content (n={len(sample_gc):,})")

    plt.tight_layout()
    plt.savefig(output_dir / "nanopore_analysis.png", dpi=300, bbox_inches="tight")
    plt.show()

    print(f"📊 Visualizations saved to {output_dir / 'nanopore_analysis.png'}")


def print_summary(stats, composition, homopolymer_data):
    """
    Print a comprehensive summary of the analysis.
    """

    print("\n" + "=" * 60)
    print("🧬 NANOPORE SEQUENCING DATA ANALYSIS SUMMARY")
    print("=" * 60)

    print("\n📊 BASIC STATISTICS:")
    print(f"  Total reads:        {stats['total_reads']:,}")
    print(f"  Total bases:        {stats['total_bases']:,} ({stats['total_bases'] / 1e6:.1f}M)")
    print(f"  Mean read length:   {stats['mean_length']:,.0f} bp")
    print(f"  Median read length: {stats['median_length']:,.0f} bp")
    print(f"  N50:                {stats['n50']:,.0f} bp")
    print(f"  Longest read:       {stats['max_length']:,} bp")
    print(f"  Shortest read:      {stats['min_length']:,} bp")

    print("\n🧮 NUCLEOTIDE COMPOSITION:")
    for base in ["A", "T", "G", "C"]:
        pct = composition["composition_percent"][base]
        print(f"  {base}: {pct:.1f}%")

    gc_pct = composition["composition_percent"]["G"] + composition["composition_percent"]["C"]
    print(f"  GC content: {gc_pct:.1f}% (mean per read: {composition['mean_gc_content']:.1f}%)")

    print("\n🔄 HOMOPOLYMER ANALYSIS:")
    total_homopolymers = sum(homopolymer_data[base]["count"] for base in ["A", "T", "G", "C"])
    print(f"  Total homopolymers ≥4bp: {total_homopolymers:,}")

    for base in ["A", "T", "G", "C"]:
        data = homopolymer_data[base]
        if data["count"] > 0:
            print(
                f"  {base}: {data['count']:,} runs (max: {data['max_length']}bp, avg: {data['mean_length']:.1f}bp)"
            )

    print("\n💡 ML ENGINEERING INSIGHTS:")
    print(f"  - Variable sequence lengths: {stats['min_length']}-{stats['max_length']} bp")
    print("  - This is like having variable-length text documents")
    print(f"  - GC bias: {abs(50 - gc_pct):.1f}% deviation from expected 50%")
    print("  - Homopolymer errors are the main ML challenge for nanopore data")
    print("  - Long reads enable genome assembly but have higher error rates")


def main():
    """
    Main analysis pipeline.
    """

    # Check if data exists
    data_file = Path("data/raw/ecoli_nanopore_reads.fastq.gz")
    if not data_file.exists():
        print("❌ Data not found! Run download_data.py first:")
        print("   uv run python scripts/download_data.py")
        return

    print("🧬 Starting Nanopore Read Analysis")
    print("=" * 50)

    # Load reads
    reads_data = load_reads(data_file)

    # Basic statistics
    print("\n📊 Calculating basic statistics...")
    stats = calculate_basic_stats(reads_data)

    # Sequence composition
    composition = analyze_sequence_composition(reads_data)

    # Homopolymer analysis
    homopolymer_data = analyze_homopolymers(reads_data)

    # Create visualizations
    print("\n📈 Creating visualizations...")
    create_visualizations(stats, composition, homopolymer_data, reads_data)

    # Print summary
    print_summary(stats, composition, homopolymer_data)

    # Save data for notebook analysis
    output_dir = Path("results")

    # Save summary statistics
    summary_data = {
        "basic_stats": stats,
        "composition": composition,
        "homopolymers": homopolymer_data,
    }

    import json

    with open(output_dir / "analysis_summary.json", "w") as f:
        # Convert numpy types to native Python types for JSON serialization
        def convert_numpy(obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, np.number):
                return obj.item()
            elif isinstance(obj, dict):
                return {key: convert_numpy(value) for key, value in obj.items()}
            elif isinstance(obj, list):
                return [convert_numpy(item) for item in obj]
            else:
                return obj

        json.dump(convert_numpy(summary_data), f, indent=2)

    print(f"\n💾 Analysis data saved to {output_dir / 'analysis_summary.json'}")

    print("""
🎉 Analysis Complete!

Next steps:
1. Open the Jupyter notebook: uv run jupyter lab
2. Explore the visualizations in results/nanopore_analysis.png
3. Use the analysis_summary.json for further ML work

Happy genomics! 🧬
    """)


if __name__ == "__main__":
    main()
