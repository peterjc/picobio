#!/usr/bin/env python
import os
import sys

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
    for fasta in fasta_list:
        ispcr_file = os.path.join(tmp, os.path.basename(fasta) + ".tsv")
        if not os.path.isfile(ispcr_file):
            print(f"Calling isPrc on {fasta}")
            cmd = f"isPcr '{fasta}' '{primer_file}' '{ispcr_file}'"
            if os.system(cmd):
                sys.exit("ERROR: Calling ispcr failed\n")


if __name__ == "__main__":
    main()
