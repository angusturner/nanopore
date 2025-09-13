#!/usr/bin/env python3
"""
Verify downloaded nanopore sequencing data.

This script checks that all required data files exist, are the correct size,
and can be properly opened and parsed.
"""

import gzip
import sys
from pathlib import Path

from Bio import SeqIO


def check_file_exists(file_path, description, expected_min_size_mb=0):
    """Check if file exists and meets minimum size requirements"""
    path = Path(file_path)

    if not path.exists():
        print(f"❌ {description}: Missing!")
        print(f"   Expected location: {file_path}")
        return False

    # Check file size
    size_mb = path.stat().st_size / (1024 * 1024)
    print(f"✅ {description}: {path.name} ({size_mb:.1f} MB)")

    if expected_min_size_mb > 0 and size_mb < expected_min_size_mb:
        print(f"⚠️  Warning: File seems small (expected >{expected_min_size_mb}MB)")
        return False

    return True


def verify_fastq_file(file_path):
    """Verify FASTQ file can be opened and parsed"""
    print("\n🔍 Verifying FASTQ file structure...")

    try:
        # Handle compressed files
        open_func = gzip.open if str(file_path).endswith(".gz") else open
        open_mode = "rt" if str(file_path).endswith(".gz") else "r"

        with open_func(file_path, open_mode) as handle:
            # Try to read first few records
            record_count = 0
            total_bases = 0

            for record in SeqIO.parse(handle, "fastq"):
                record_count += 1
                total_bases += len(record.seq)

                # Check first record has required components
                if record_count == 1:
                    if not record.id:
                        print("❌ First record missing ID")
                        return False
                    if not record.seq:
                        print("❌ First record missing sequence")
                        return False
                    if (
                        not hasattr(record, "letter_annotations")
                        or "phred_quality" not in record.letter_annotations
                    ):
                        print("❌ First record missing quality scores")
                        return False

                    print(f"   Sample read ID: {record.id[:50]}...")
                    print(f"   Sample sequence: {str(record.seq)[:50]}...")
                    print(
                        f"   Quality scores: {record.letter_annotations['phred_quality'][:10]}..."
                    )

                # Stop after checking 100 records for efficiency
                if record_count >= 100:
                    break

        if record_count == 0:
            print("❌ No valid FASTQ records found")
            return False

        avg_read_length = total_bases / record_count if record_count > 0 else 0
        print("✅ FASTQ file structure valid")
        print(f"   Checked {record_count} records")
        print(f"   Average read length: {avg_read_length:.0f} bp")
        print(f"   Total bases checked: {total_bases:,}")

        return True

    except Exception as e:
        print(f"❌ Error reading FASTQ file: {e}")
        return False


def verify_fasta_file(file_path):
    """Verify FASTA file can be opened and parsed"""
    print("\n🔍 Verifying FASTA file structure...")

    try:
        with open(file_path) as handle:
            record_count = 0
            total_bases = 0

            for record in SeqIO.parse(handle, "fasta"):
                record_count += 1
                total_bases += len(record.seq)

                # Check first record
                if record_count == 1:
                    if not record.id:
                        print("❌ First record missing ID")
                        return False
                    if not record.seq:
                        print("❌ First record missing sequence")
                        return False

                    print(f"   Genome ID: {record.id}")
                    print(f"   Sequence preview: {str(record.seq)[:50]}...")

                # Most reference genomes are single sequence files
                if record_count > 10:
                    break

        if record_count == 0:
            print("❌ No valid FASTA records found")
            return False

        print("✅ FASTA file structure valid")
        print(f"   Contains {record_count} sequence(s)")
        print(f"   Total bases: {total_bases:,}")

        return True

    except Exception as e:
        print(f"❌ Error reading FASTA file: {e}")
        return False


def check_data_compatibility():
    """Check that the data files are compatible with analysis scripts"""
    print("\n🔗 Checking data compatibility...")

    # Check that analysis script exists and mentions the correct filenames
    analyze_script = Path("scripts/analyze_reads.py")
    if analyze_script.exists():
        print("✅ Analysis script found")
    else:
        print("⚠️  Analysis script not found")
        return False

    # Check notebook exists
    notebook = Path("notebooks/nanopore_exploration.ipynb")
    if notebook.exists():
        print("✅ Exploration notebook found")
    else:
        print("⚠️  Exploration notebook not found")
        return False

    print("✅ All analysis tools available")
    return True


def main():
    """Main verification workflow"""
    print("🔍 Nanopore Data Verification Script")
    print("=" * 50)

    all_good = True

    # Files to check with expected minimum sizes
    files_to_check = [
        ("data/raw/ecoli_nanopore_reads.fastq.gz", "Nanopore reads (FASTQ)", 100),  # Expect >100MB
        ("data/references/ecoli_reference.fasta", "Reference genome (FASTA)", 4),  # Expect >4MB
    ]

    print("📋 Checking file existence and sizes...")

    # Basic file checks
    for file_path, description, min_size in files_to_check:
        if not check_file_exists(file_path, description, min_size):
            all_good = False

    # If basic checks pass, do deeper validation
    if all_good:
        print("\n📖 Performing detailed file validation...")

        # Verify FASTQ file
        fastq_path = Path("data/raw/ecoli_nanopore_reads.fastq.gz")
        if fastq_path.exists() and not verify_fastq_file(fastq_path):
            all_good = False

        # Verify FASTA file
        fasta_path = Path("data/references/ecoli_reference.fasta")
        if fasta_path.exists() and not verify_fasta_file(fasta_path):
            all_good = False

        # Check compatibility
        if not check_data_compatibility():
            all_good = False

    # Final report
    print("\n" + "=" * 50)
    if all_good:
        print("🎉 All data verification checks passed!")
        print("""
✅ Your data is ready for analysis!

Next steps:
1. Run basic analysis:
   uv run python scripts/analyze_reads.py

2. Start interactive exploration:
   uv run jupyter lab notebooks/nanopore_exploration.ipynb

3. Results will be saved to results/ directory

💡 Data Summary:
- FASTQ file contains nanopore reads with quality scores
- FASTA file contains the reference genome
- Both files are properly formatted and readable
- Analysis scripts are ready to use
        """)
    else:
        print("❌ Data verification failed!")
        print("""
Issues detected. Please:

1. Re-run the download script:
   uv run python scripts/download_data.py

2. Check your internet connection

3. Ensure you have sufficient disk space

4. If problems persist, check the file URLs in the download script
        """)
        sys.exit(1)


if __name__ == "__main__":
    main()
