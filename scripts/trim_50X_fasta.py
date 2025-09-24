#!/usr/bin/env python3
"""
trim_fasta.py

Trim a fixed number of bases from the 5' and 3' ends of nucleotide sequences in a FASTA file.

This script is useful for removing sequencing adapters, primers, or padding regions from synthetic constructs.

The intended use is to remove 50X (e.g., 504, 506, 507, etc.) primers from reference FASTA files

Usage:
    python scripts/trim_50X_fasta.py input.fa output_trimmed.fa
    python scripts/trim_50X_fasta.py input.fa output_trimmed.fa --trim5 25 --trim3 15

    
Example:
    python scripts/trim_50X_fasta.py refs/Twist_384-1536_Unique/504_Primers/unique_overlap_lib384_gene_full_with_504.fasta refs/Twist_384-1536_Unique/unique_overlap_lib384_gene_full_wo_504.fasta
    python scripts/trim_50X_fasta.py refs/Twist_384-1536_Unique/504_Primers/unique_overlap_lib1536_gene_full_with_504.fasta refs/Twist_384-1536_Unique/unique_overlap_lib1536_gene_full_wo_504.fasta

    python scripts/trim_50X_fasta.py refs/AlignFdn/504_Primers/AF.lib2.384_wPad_wRE_wPrim.genes refs/AlignFdn/AF.lib2.384_wPad_wRE_n504.fasta
    python scripts/trim_50X_fasta.py refs/AlignFdn/504_Primers/AF.lib16.1536_wPad_wRE_wPrim.genes refs/AlignFdn/AF.lib16.1536_wPad_wRE_n504.fasta
    python scripts/trim_50X_fasta.py refs/AlignFdn/504_Primers/AF.lib18.1536_wPad_wRE_wPrim.genes refs/AlignFdn/AF.lib18.1536_wPad_wRE_n504.fasta
    python scripts/trim_50X_fasta.py refs/AlignFdn/504_Primers/AF.lib20.1536_wPad_wRE_wPrim.genes refs/AlignFdn/AF.lib20.1536_wPad_wRE_n504.fasta

Arguments:
    input.fa            Input FASTA file containing raw sequences.
    output_trimmed.fa   Output FASTA file to write the trimmed sequences.
    --trim5             Number of bases to trim from the 5' end (default: 20).
    --trim3             Number of bases to trim from the 3' end (default: 20).

Requirements:
    - Python 3
    - Biopython (`pip install biopython`)
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse

def trim_fasta(input_fasta, output_fasta, trim5=20, trim3=20):
    trimmed_records = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        trimmed_seq = record.seq[trim5:len(record.seq) - trim3]
        trimmed_record = SeqRecord(
            trimmed_seq,
            id=record.id,
            description=record.description
        )
        trimmed_records.append(trimmed_record)

    SeqIO.write(trimmed_records, output_fasta, "fasta")
    print(f"Trimmed {len(trimmed_records)} sequences and saved to {output_fasta}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Trim 5' and 3' ends from FASTA sequences.")
    parser.add_argument("input_fasta", help="Path to input FASTA file")
    parser.add_argument("output_fasta", help="Path to output trimmed FASTA file")
    parser.add_argument("--trim5", type=int, default=20, help="Number of bases to trim from 5' end (default: 20)")
    parser.add_argument("--trim3", type=int, default=20, help="Number of bases to trim from 3' end (default: 20)")
    args = parser.parse_args()

    trim_fasta(args.input_fasta, args.output_fasta, args.trim5, args.trim3)
