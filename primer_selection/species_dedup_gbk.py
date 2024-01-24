#!/usr/bin/env python
r"""De-duplicate GenBank file to species level with FASTA output.

Parses stdin, writes to stdout. Takes the first entry of each
species found (discarding duplicates). The output FASTA files
are named ">{accession} {lineage}\n{sequence}\n" with the
taxonomic lineage from the GenBank header.
"""
import sys

from Bio import SeqIO
from Bio.Seq import UndefinedSequenceError

ignored = 0
unique = 0
counts = {}
for record in SeqIO.parse(sys.stdin, "genbank"):
    lineage = (
        ";".join(record.annotations["taxonomy"]) + ";" + record.annotations["organism"]
    )
    try:
        counts[lineage] += 1
    except KeyError:
        try:
            sys.stdout.write(f">{record.id} {lineage}\n{record.seq}\n")
            unique += 1
            counts[lineage] = 1
        except UndefinedSequenceError:
            sys.stderr.write(f"WARNING - Ignoring {record.id} as no sequence\n")

sys.stderr.write(
    f"{unique} unique counts in {sum(counts.values())} sequences; plus {ignored} ignored\n"
)
# for sp, count in sorted(counts.items()):
#    sys.stderr.write(f"{sp}\t{count}\n")
