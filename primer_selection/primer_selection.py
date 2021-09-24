#!/usr/bin/env python
import os
import sys

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


def load_isprc(isprc_filename, primer_hits):
    """Parse the FASTA output from Jim Kent's isPcr tool.

    Adds (reference, sequence) entries for (left, right) in primer_hits dict.
    """
    ref_name = os.path.splitext(os.path.split(isprc_filename)[1])[0]
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
    if len(sys.argv) < 3:
        sys.exit(
            "ERROR: At least two arguments required, primer TSV and one or more FASTA"
        )

    primers = load_primers(sys.argv[1])
    if not primers:
        sys.exit(f"ERROR: No primers identified in {sys.argv[1]}")
    primer_file = os.path.join(tmp, "primers.tsv")
    with open(primer_file, "w") as handle:
        for name, left, right in primers:
            handle.write(f"{name}\t{left}\t{right}\n")

    fasta_list = sys.argv[2:]
    primer_hits = {(left, right): [] for name, left, right in primers}
    for fasta in fasta_list:
        ispcr_file = os.path.join(tmp, os.path.basename(fasta) + ".tsv")
        if not os.path.isfile(ispcr_file):
            print(f"Calling isPrc on {fasta}")
            cmd = f"isPcr '{fasta}' '{primer_file}' '{ispcr_file}'"
            if os.system(cmd):
                sys.exit("ERROR: Calling ispcr failed\n")
        load_isprc(ispcr_file, primer_hits)

    for left, right in primer_hits:
        print(left, right, primer_hits[left, right])


if __name__ == "__main__":
    main()
