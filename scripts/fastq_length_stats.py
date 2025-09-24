#!/usr/bin/env python3
"""
Calculate mean, median, min, and max read lengths from a FASTQ file (.fastq or .fastq.gz).

Also saves the output to a .log file in the logs/ directory.

Usage:
    python fastq_length_stats.py input.fastq[.gz]
"""

import sys
import gzip
import statistics
import os

def fastq_lengths(fastq_file):
    """Generator that yields sequence lengths from a FASTQ file."""
    opener = gzip.open if fastq_file.endswith(".gz") else open
    mode = "rt" if fastq_file.endswith(".gz") else "r"

    with opener(fastq_file, mode) as f:
        for i, line in enumerate(f):
            if i % 4 == 1:  # sequence line (2nd of every 4 lines)
                yield len(line.strip())

def main():
    if len(sys.argv) != 2:
        print("Usage: python fastq_length_stats.py input.fastq[.gz]")
        sys.exit(1)

    fastq_file = sys.argv[1]
    lengths = list(fastq_lengths(fastq_file))

    if not lengths:
        print("No sequences found in the FASTQ file.")
        sys.exit(1)

    mean_len = statistics.mean(lengths)
    median_len = statistics.median(lengths)
    min_len = min(lengths)
    max_len = max(lengths)

    # Prepare output text
    output_lines = [
        f"File: {fastq_file}",
        f"Number of reads: {len(lengths)}",
        f"Mean length:   {mean_len:.2f}",
        f"Median length: {median_len:.2f}",
        f"Min length:    {min_len}",
        f"Max length:    {max_len}",
    ]
    output_text = "\n".join(output_lines)

    # Print to terminal
    print(output_text)

if __name__ == "__main__":
    main()
