#!/usr/bin/env python3

import argparse
import os
from collections import Counter
from Bio import SeqIO

def truncate_at_whitespace(header: str) -> str:
    """
    Return everything up to the first whitespace.
    If the header starts with '@', remove that first.
    """
    if header.startswith('@'):
        header = header[1:]
    return header.split()[0]

def parse_sam(sam_path):
    ref_counts = Counter()
    read_to_ref = {}
    with open(sam_path, 'r') as sam_in:
        for line in sam_in:
            if line.startswith('@'):
                continue
            f = line.split('\t')
            raw_id = f[0]; flag = int(f[1]); ref = f[2]
            if ref == '*' or (flag & 256) or (flag & 2048):
                continue
            rid = truncate_at_whitespace(raw_id)
            if rid not in read_to_ref:
                read_to_ref[rid] = ref
                ref_counts[ref] += 1
    return ref_counts, read_to_ref

def write_counts(ref_counts, csv_path):
    with open(csv_path, 'w') as out:
        out.write("reference,read_count\n")
        for r,c in ref_counts.most_common():
            out.write(f"{r},{c}\n")

def split_fasta(fasta_in, read_to_ref, top_refs, outdir):
    os.makedirs(outdir, exist_ok=True)
    handles = {r: open(os.path.join(outdir, f"{r}_reads.fasta"), 'w')
               for r in top_refs}
    written = Counter()
    for rec in SeqIO.parse(fasta_in, "fasta"):
        rid = rec.id
        if rid in read_to_ref and read_to_ref[rid] in handles:
            SeqIO.write(rec, handles[read_to_ref[rid]], "fasta")
            written[read_to_ref[rid]] += 1
    for h in handles.values():
        h.close()
    return written

def main():
    p = argparse.ArgumentParser(
        description="Split consensus FASTA by top-N SAM alignments")
    p.add_argument("-s","--sam",      required=True, help="Input SAM")
    p.add_argument("-f","--fasta",    required=True, help="Consensus FASTA")
    p.add_argument("-c","--counts",   required=True, help="Output CSV")
    p.add_argument("-o","--outdir",   required=True, help="Output dir for per-ref FASTAs")
    p.add_argument("-t","--top", type=int, default=4, help="Number of top refs")
    args = p.parse_args()

    # Step 1
    ref_counts, read_to_ref = parse_sam(args.sam)
    total = sum(ref_counts.values())
    print(f"Total primary mapped reads: {total}")

    # Step 2
    write_counts(ref_counts, args.counts)
    print(f"Wrote counts → {args.counts}")

    # Step 3
    top_refs = [r for r,_ in ref_counts.most_common(args.top)]
    print(f"Top {args.top} refs:")
    for r in top_refs:
        print(f"  {r} ({ref_counts[r]} reads)")

    # Step 4
    written = split_fasta(args.fasta, read_to_ref, top_refs, args.outdir)
    print(f"Wrote per-ref FASTAs into {args.outdir}:")
    for r in top_refs:
        print(f"  {r}_reads.fasta: {written[r]} seqs")

if __name__=="__main__":
    main()
