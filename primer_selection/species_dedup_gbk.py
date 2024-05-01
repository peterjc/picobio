#!/usr/bin/env python
r"""De-duplicate GenBank file to species level with FASTA output.

Parses stdin, writes to stdout. Takes the first entry of each
species found (discarding duplicates). The output FASTA files
are named ">{accession} {lineage}\n{sequence}\n" with the
taxonomic lineage from the GenBank header.
"""
import argparse
import sys

from Bio import SeqIO
from Bio.Seq import UndefinedSequenceError


parser = argparse.ArgumentParser(
    prog="species_dedup_gbk.py",
    description=("Depulicate GenBank files to species level FASTA file(s)."),
)
parser.add_argument(
    "filenames",
    type=argparse.FileType("r"),
    nargs="*",
    metavar="GBK",
    default=[sys.stdin],
    help="One or more GenBank files (default stdin)",
)
parser.add_argument(
    "-o",
    "--output",
    type=str,
    metavar="STEM",
    help="Output filename stem for generating reports, default all to stdout.",
)
parser.add_argument(
    "-c",
    "--chunk",
    type=int,
    metavar="COUNT",
    default=10000,
    help="Number of sequences per file if using the output stemp option.",
)
options = parser.parse_args()

if options.output:
    filename_pattern = options.output + ".0000.fasta"
    out_handle = open(filename_pattern, "w")
else:
    filename_pattern = None
    out_handle = sys.stdout

# if not options.filenames:
#    sys.stderr.write("Parsing from stdin...\n")
#    options.filenames = [sys.stdin]

ignored = 0
unique = 0
counts = {}
chunk = 0
for in_handle in options.filenames:
    sys.stderr.write(f"Parsing {in_handle.name}\n")
    for record in SeqIO.parse(in_handle, "genbank"):
        lineage = (
            ";".join(record.annotations["taxonomy"])
            + ";"
            + record.annotations["organism"]
        )
        try:
            counts[lineage] += 1
        except KeyError:
            try:
                out_handle.write(f">{record.id} {lineage}\n{record.seq}\n")
                unique += 1
                counts[lineage] = 1
                if filename_pattern and unique % options.chunk == 0:
                    out_handle.close()
                    sys.stderr.write(f"{unique} so far, finished {out_handle.name}\n")
                    chunk += 1
                    filename_pattern = f"{options.output}.{chunk:04d}.fasta"
                    out_handle = open(filename_pattern, "w")
            except UndefinedSequenceError:
                sys.stderr.write(f"WARNING - Ignoring {record.id} as no sequence\n")

if filename_pattern:
    out_handle.close()
    sys.stderr.write(f"{unique} unique, finished {out_handle.name}\n")

sys.stderr.write(
    f"{unique} unique counts in {sum(counts.values())} sequences; plus {ignored} ignored\n"
)
# for sp, count in sorted(counts.items()):
#    sys.stderr.write(f"{sp}\t{count}\n")
