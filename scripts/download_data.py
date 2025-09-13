#!/usr/bin/env python3
"""
Download real nanopore sequencing data for analysis.

This script downloads publicly available E. coli nanopore reads and reference genome.
Perfect for learning bioinformatics without needing actual sequencing hardware!
"""

import sys
import time
from pathlib import Path

import requests


def download_file(url, output_path, description=""):
    """Download a file with resume capability and automatic retries"""
    print(f"Downloading {description}...")
    print(f"URL: {url}")
    print(f"Output: {output_path}")

    # Strategy 1: Resume downloads - check existing file
    resume_byte_pos = 0
    if output_path.exists():
        resume_byte_pos = output_path.stat().st_size
        print(f"📄 Found partial file, resuming from byte {resume_byte_pos:,}")

    # Strategy 2: Automatic retries with exponential backoff
    max_retries = 5
    for attempt in range(max_retries):
        try:
            headers = {}
            if resume_byte_pos > 0:
                headers["Range"] = f"bytes={resume_byte_pos}-"

            response = requests.get(url, stream=True, headers=headers, timeout=30)
            response.raise_for_status()

            # Get total file size (accounting for resume)
            if "content-range" in response.headers:
                # Partial content response
                total_size = int(response.headers["content-range"].split("/")[-1])
            else:
                # Full file download
                total_size = int(response.headers.get("content-length", 0))

            downloaded_size = resume_byte_pos

            # Open in append mode if resuming, write mode otherwise
            mode = "ab" if resume_byte_pos > 0 else "wb"

            with open(output_path, mode) as f:
                for chunk in response.iter_content(chunk_size=32768):  # Larger chunks = faster
                    if chunk:
                        f.write(chunk)
                        downloaded_size += len(chunk)

                        if total_size > 0:
                            progress = (downloaded_size / total_size) * 100
                            print(
                                f"\rProgress: {progress:.1f}% ({downloaded_size:,} / {total_size:,} bytes) "
                                f"[Attempt {attempt + 1}/{max_retries}]",
                                end="",
                            )

            print(f"\n✅ Downloaded {output_path.name} ({downloaded_size:,} bytes)")
            return True

        except (requests.RequestException, OSError) as e:
            print(f"\n⚠️  Download failed (attempt {attempt + 1}/{max_retries}): {e}")

            if attempt < max_retries - 1:
                # Strategy 2: Exponential backoff - wait longer each time
                wait_time = 2**attempt  # 1, 2, 4, 8 seconds
                print(f"Retrying in {wait_time} seconds...")
                time.sleep(wait_time)

                # Update resume position in case we got some data
                if output_path.exists():
                    resume_byte_pos = output_path.stat().st_size
            else:
                print(f"❌ All {max_retries} attempts failed for {description}")
                return False

        except Exception as e:
            print(f"\n❌ Unexpected error downloading {description}: {e}")
            return False

    return False


def setup_directories():
    """Create necessary directories"""
    dirs = [Path("data/raw"), Path("data/references"), Path("results")]

    for dir_path in dirs:
        dir_path.mkdir(parents=True, exist_ok=True)
        print(f"📁 Created directory: {dir_path}")


def download_ecoli_data():
    """
    Download E. coli nanopore sequencing data.

    Why E. coli?
    - Small genome (~4.6M bases) - manageable for learning
    - Well-characterized organism (K-12 MG1655 reference strain)
    - Common model organism in bioinformatics
    - This dataset is from Zenodo with real MinION data
    """

    # E. coli nanopore reads from Zenodo - basecalled FASTQ
    reads_url = (
        "https://zenodo.org/records/7995806/files/bonito_local_basecalled.fastq.gz?download=1"
    )
    reads_path = Path("data/raw/ecoli_nanopore_reads.fastq.gz")

    success = download_file(
        reads_url, reads_path, "E. coli nanopore reads (FASTQ format with quality scores)"
    )

    if success:
        print("""
📊 About this data:
- E. coli strain K-12 MG1655 (reference strain)
- Oxford Nanopore MinION sequencer (Flow Cell R9.4.1)
- Basecalled with Bonito (high-quality basecaller)
- FASTQ format (sequences + quality scores!)
- Real nanopore data with authentic error patterns
- Perfect for ML engineers learning bioinformatics
        """)

    return success


def download_reference_genome():
    """
    Download E. coli reference genome from the same Zenodo dataset.

    What's a reference genome?
    - The "correct" sequence we can compare our reads against
    - Used for alignment, variant calling, assembly validation
    - This is E. coli strain K-12 MG1655 - the most studied E. coli strain
    - FROM THE SAME DATASET as our sequencing reads!
    """

    ref_url = "https://zenodo.org/records/7995806/files/Ecoli_K12_MG1655.fasta?download=1"
    ref_path = Path("data/references/ecoli_reference.fasta")

    success = download_file(ref_url, ref_path, "E. coli reference genome (K-12 MG1655)")

    if success:
        print(f"✅ Reference genome saved as: {ref_path}")
        print("""
🧬 About the reference:
- E. coli K-12 MG1655 complete genome
- Same strain as the sequencing data!
- From the same Zenodo dataset
- Perfect match for alignment and analysis
- 4,641,652 base pairs, single circular chromosome
        """)

    return success


def main():
    """Main download workflow"""
    print("🧬 Nanopore Data Download Script")
    print("=" * 50)

    # Setup
    setup_directories()

    # Download data
    print("\n📥 Downloading data...")
    ecoli_success = download_ecoli_data()
    ref_success = download_reference_genome()

    if ecoli_success and ref_success:
        print("""
🎉 Download complete!

Next steps:
1. Verify your downloads:
   uv run python scripts/verify_data.py

2. Run basic analysis:
   uv run python scripts/analyze_reads.py

3. Start interactive notebook:
   uv run jupyter lab notebooks/

💡 ML Engineer Tips:
- These are real sequencing reads with authentic error patterns
- Perfect for learning about sequence data before building ML models
- Each read represents a single DNA molecule that was sequenced
        """)
    else:
        print("\n❌ Downloads failed. Check your internet connection and try again.")
        sys.exit(1)


if __name__ == "__main__":
    main()
