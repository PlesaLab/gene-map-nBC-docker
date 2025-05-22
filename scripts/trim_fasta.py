#!/usr/bin/env python3

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
