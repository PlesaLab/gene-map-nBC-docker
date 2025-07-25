#!/usr/bin/env python3
"""
extract_seqs.py

Extracts sequence fragments between two specified restriction enzyme 
(RE) sites (including the RE sites themselves) from all records in 
an input FASTA file. Also searches for the reverse complement 
orientation of the RE sites and extracts those fragments 
(reverse-complemented back to forward orientation).

Usage:
    python extract_seqs.py input.fasta output.fasta

Outputs:
  - output.fasta: FASTA file containing extracted fragments.
  - output.log: Log file recording counts of forward and reverse-complement matches.
"""
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input.fasta output.fasta")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    log_file = output_file.rsplit('.', 1)[0] + '.log'

    # Define RE site sequences (CURRENT: AlignFdn)
    start_site = "GCATGGTCTCATCGGCACCTGCATGAATCT"
    end_site   = "TAATTACTGCAGGTGGACCTGAGACCGCAT"
    # Compute reverse complements
    rc_start = str(Seq(start_site).reverse_complement())
    rc_end   = str(Seq(end_site).reverse_complement())

    forward_count = 0
    reverse_count = 0
    extracted = []

    # Parse input FASTA
    for record in SeqIO.parse(input_file, "fasta"):
        seq_str = str(record.seq)
        # Forward orientation
        f_start = seq_str.find(start_site)
        if f_start != -1:
            f_end = seq_str.find(end_site, f_start + len(start_site))
            if f_end != -1:
                fragment = seq_str[f_start : f_end + len(end_site)]
                new_rec = SeqRecord(Seq(fragment), id=record.id, description="")
                extracted.append(new_rec)
                forward_count += 1
                continue
        # Reverse orientation
        r_start = seq_str.find(rc_start)
        if r_start != -1:
            r_end = seq_str.find(rc_end, r_start + len(rc_start))
            if r_end != -1:
                fragment = seq_str[r_start : r_end + len(rc_end)]
                # Reverse-complement back to forward orientation
                fragment_rc = str(Seq(fragment).reverse_complement())
                new_rec = SeqRecord(Seq(fragment_rc), id=record.id + "_rc", description="")
                extracted.append(new_rec)
                reverse_count += 1

    # Write extracted sequences
    SeqIO.write(extracted, output_file, "fasta")

    # Write log file
    with open(log_file, 'w') as log:
        log.write(f"Forward matches: {forward_count}\n")
        log.write(f"Reverse-complement matches: {reverse_count}\n")

    print(f"Done: {forward_count} forward, {reverse_count} reverse matches.\n"
          f"Results: {output_file}\nLog: {log_file}")

if __name__ == '__main__':
    main()
