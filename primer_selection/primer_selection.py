#!/usr/bin/env python
import argparse
import os
import sys
from collections import defaultdict
from string import ascii_uppercase

from Bio.SeqIO.FastaIO import SimpleFastaParser

tmp = "/tmp/primer_selection"
if not os.path.isdir(tmp):
    os.mkdir(tmp)


def load_primers(tsv_filename):
    """Load primer TSV file into list of 3-tuples."""
    answer = []
    with open(tsv_filename) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) == 2:
                left, right = parts
                name = f"P{len(answer)}"
            else:
                name, left, right = parts[:3]
            answer.append((name, left, right))
    return answer


def load_isprc(isprc_filename, ref_name, primer_hits):
    """Parse the FASTA output from Jim Kent's isPcr tool.

    Adds (reference, sequence) entries for (left, right) in primer_hits dict.
    """
    with open(isprc_filename) as handle:
        for title, seq in SimpleFastaParser(handle):
            amp_region, name, length, left, right = title.split()
            assert (
                left,
                right,
            ) in primer_hits, f"Stale cache? Why {name}, {left}, {right}"
            assert (
                length == f"{len(seq)}bp"
            ), f"Expected length {len(seq)} from sequence, yet {length}"
            seq = seq.upper()
            # chrom, rest = region.rsplit(":", 1)
            # if "+" in rest:
            #    strand = "+"
            #    start, end = rest.split("+")
            # else:
            #    strand = "-"
            #    start, end = rest.split("-")
            primer_hits[left, right].append((ref_name, seq))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("primers", metavar="TSV", help="TSV file of primers")
    parser.add_argument(
        "references", nargs="+", metavar="FASTA", help="One or more FASTA files"
    )
    args = parser.parse_args()

    if not args.primers:
        sys.exit("ERROR: Missing primer TSV file")
    if not args.references:
        sys.exit("ERROR: Missing FASTA file(s)")

    primers = load_primers(args.primers)
    if not primers:
        sys.exit(f"ERROR: No primers identified in {args.primers}")
    primer_file = os.path.join(tmp, "primers.tsv")
    with open(primer_file, "w") as handle:
        for name, left, right in primers:
            handle.write(f"{name}\t{left}\t{right}\n")

    fasta_list = args.references
    ref_list = [os.path.splitext(os.path.split(_)[1])[0] for _ in fasta_list]
    primer_hits = {(left, right): [] for name, left, right in primers}
    for fasta, ref_name in zip(fasta_list, ref_list):
        ispcr_file = os.path.join(tmp, ref_name + ".tsv")
        if not os.path.isfile(ispcr_file):
            print(f"Calling isPrc on {fasta}")
            cmd = f"isPcr '{fasta}' '{primer_file}' '{ispcr_file}'"
            if os.system(cmd):
                sys.exit("ERROR: Calling ispcr failed\n")
        load_isprc(ispcr_file, ref_name, primer_hits)

    # print(f"Have {len(fasta_list)} references")
    amplicons = defaultdict(dict)
    for (left, right), values in primer_hits.items():
        for ref_name, seq in values:
            try:
                amplicons[left, right][seq] += 1
            except KeyError:
                amplicons[left, right][seq] = 1

    # Assign letters to each unique sequnce for each amplicon: A, B, ...
    amplicon_alias = {}
    for (left, right), seq_counts in amplicons.items():
        for i, (count, seq) in enumerate(
            reversed(sorted((count, seq) for seq, count in seq_counts.items()))
        ):
            amplicon_alias[left, right, seq] = ascii_uppercase[i]

    print(f"Have {len(ref_list)} references, and {len(primer_hits)} primers")
    for ref in ref_list:
        print(
            ";".join(
                "".join(
                    amplicon_alias[left, right, s]
                    for (r, s) in primer_hits[left, right]
                    if r == ref
                )
                for left, right in primer_hits
            ),
            ref,
        )


if __name__ == "__main__":
    main()
