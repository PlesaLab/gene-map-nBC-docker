#!/usr/bin/env python3
import os
import regex  # for approximate matching
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import collections
import csv
import random
import argparse

def DNA(length):
    return ''.join(random.choice('CGTA') for _ in range(length))

def open_maybe_gz(filename):
    return gzip.open(filename, 'rt') if filename.endswith('.gz') else open(filename, 'r')

def main():
    parser = argparse.ArgumentParser(
        description="Process barcodes from a FASTQ file without requiring a config file."
    )
    parser.add_argument(
        'input_file',
        help='Path to the input FASTQ or FASTQ.GZ file'
    )
    parser.add_argument(
        'output_prefix',
        help='Output prefix (directory + base name) for generated files'
    )
    parser.add_argument(
        '--startsitelength',
        type=int,
        required=True,
        help='Number of bases in the start_site (for slicing)'
    )
    parser.add_argument(
        '--barcode_length',
        type=int,
        required=True,
        help='Length (in bp) of the random barcode to generate'
    )
    parser.add_argument(
        '--motif',
        required=True,
        help='Regex for flanking constant sequences around the barcode (e.g. ((TACGAATCGGGG)(.)(TGGTAACTAACG)){e<2})'
    )
    parser.add_argument(
        '--start_site',
        required=True,
        help='5\' flanking sequence for barcode start site'
    )
    parser.add_argument(
        '--end_site',
        required=True,
        help='3\' flanking sequence for barcode end site'
    )
    args = parser.parse_args()

    input_file = args.input_file
    output_prefix = args.output_prefix
    startsitelength = args.startsitelength
    barcode_length = args.barcode_length
    motif = args.motif
    start_site = args.start_site
    end_site = args.end_site

    out_dir = os.path.dirname(output_prefix) or '.'
    logs_dir = 'logs'
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    with open_maybe_gz(input_file) as handle:
        records = list(SeqIO.parse(handle, 'fastq'))
    total = len(records)
    print(f"{total} total records")

    rc_count = nomatchcount = count_no_startsite = count_no_endsite = bafframecount = goodaa = 0
    d = collections.defaultdict(list)
    dnum = {}
    outputseqs_noN = []

    for i, record in enumerate(records, 1):
        if i % 50000 == 0:
            print(f"On record: {i}")
        seq = str(record.seq).upper()
        match = regex.search(motif, seq, regex.BESTMATCH)
        if not match:
            seq = str(record.seq.reverse_complement()).upper()
            match = regex.search(motif, seq, regex.BESTMATCH)
            rc_count += 1

        if match:
            bc = DNA(barcode_length)
            d[bc].append(seq)
            if bc not in dnum:
                dnum[bc] = 1
                si = seq.find(start_site)
                if si != -1:
                    ei = seq.find(end_site)
                    if ei != -1:
                        subseq = seq[si + startsitelength:ei]
                        if 'N' not in subseq:
                            outputseqs_noN.append(SeqRecord(id=bc, seq=Seq(subseq), description=""))
                            dnum[bc] = dnum.get(bc, 0) + 0  # keep count
                            if len(subseq) % 3 != 0:
                                bafframecount += 1
                            aa = Seq(subseq).translate()
                            if '*' not in str(aa):
                                goodaa += 1
                    else:
                        count_no_endsite += 1
                else:
                    count_no_startsite += 1
            else:
                dnum[bc] += 1
        else:
            nomatchcount += 1

    print(f"{rc_count} times tried RC out of {total}")
    print(f"{nomatchcount} reads failed completely out of {total}")
    print(f"{len(d)} BCs found out of {total}")
    print(f"{count_no_startsite} no start site out of {len(d)}")
    print(f"{count_no_endsite} no end site out of {len(d)-count_no_startsite}")
    print(f"{len(outputseqs_noN)} no ambiguous sites out of {len(d)-count_no_startsite-count_no_endsite}")
    print(f"{goodaa} no premature stop codons out of {len(outputseqs_noN)}")

    stats_csv = f"{output_prefix}.bc_stats.csv"
    with open(stats_csv, 'w') as stats_out:
        writer = csv.writer(stats_out)
        for bc, count in dnum.items():
            writer.writerow([bc, count])

    starcode_tsv = f"{output_prefix}.bc_stats_for_starcode.tsv"
    with open(starcode_tsv, 'w') as sc_out:
        writer = csv.writer(sc_out, delimiter='\t')
        for bc, count in dnum.items():
            writer.writerow([bc, count])

    list_csv = f"{output_prefix}.bc_list.csv"
    with open(list_csv, 'w') as list_out:
        writer = csv.writer(list_out)
        for bc, seqs in d.items():
            for s in seqs:
                writer.writerow([bc, s])

    fasta_out = f"{output_prefix}.fasta"
    with open(fasta_out, 'w') as fa_out:
        SeqIO.write(outputseqs_noN, fa_out, 'fasta')

    print(f"Wrote: {stats_csv}, {starcode_tsv}, {list_csv}, {fasta_out}")

if __name__ == '__main__':
    main()
