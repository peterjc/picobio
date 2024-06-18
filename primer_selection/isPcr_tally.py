#!/usr/bin/env python
#
# Based on report_len_###.py developed for an in-silico PCR
# primer evaluation project Sep 2023 to Jan 2024.
# Renamed from isPcr_lineage_tally.py to isPcr_tally.py
# when moved from v0.9.2 to v1.0.0 and preserved FASTA IDs.
import argparse
import os
import sys

from Bio.SeqIO.FastaIO import SimpleFastaParser

if "-v" in sys.argv or "--version" in sys.argv:
    print("v1.0.0")
    sys.exit(0)

usage = """\
Parses a set of FASTA files and isPcr BED output files, to
report on the expected product legnths expected to amplify
per primer and per lineage. Produces one row for each sequence
ID in the FASTA files, with columns for the sequence ID and
description, and each primer pair. Each value is a semi-colon
separated list of any amplicon lengths for each in-silico PCR
combination. This is then used with sister script
``isPcr_lineage_report.py`` to produce tables and plots.

Inputs:

* Primer definitions in 3-column TSV format used by isPcr
* Set of reference FASTA files where the description line is
  ideally for downstream analysis a semi-colon separated taxonomic
  lineage, as produced by script ``species_dedup_gbk.py`` from the
  source lines in NCBI GenBank format files.
* Bed files (TSV) from Jim Kent's isPcr tool run on those FASTA
  files (must have matching record identifiers), and primers
  (possibly with even more primers). This may drop column 5
  (score).

The reason for ignoring column 5 is to facilitate running isPcr
in combination with script ``iupac_isPcr.py`` which expands any
IUPAC ambiguity codes in degenerate primers into a set of all
possible unambiguous interpretations of the primers. This will
result in duplicate hits differing only in score.

Example usage::

    $ sort primers.tsv | uniq | ./iupac_isPcr.py > expanded.tsv
    $ isPcr refs.fasta expanded.tsv stdout -out=bed \\
      | cut -f 1-4,6 | sort | uniq > amplicons.tsv
    $ ./isPcr_tally.py -f refs.fasta \\
      -p primers.tsv -a amplicons.tsv -o tally.tsv

Here ``primers.tsv`` is a three-column input TSV file of primer
pair name, forward, and reverse sequences -- possibly ambiguous.
Intermediate file ``expanded.tsv`` is the expanded file of
unambiguous primers, ``amplicons.tsv`` is the isPcr output in BED
format (without column 5, score), sorted and deduplicated.
Finally ``tally.tsv`` is the output tally TSV filename.
"""

parser = argparse.ArgumentParser(
    prog="isPcr_tally.py",
    description="Produce tally of Jim Kent's isPcr results, sequence vs primer.",
    epilog=usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument(
    "-f",
    "--fasta",
    nargs="+",
    metavar="FASTA",
    required=True,
    help="One or more reference FASTA files with taxonomic lineage as description.",
)
parser.add_argument(
    "-p",
    "--primers",
    nargs="+",
    metavar="TSV",
    required=True,
    help=(
        "Primer as isPcr style 3-column plain text TSV file(s) "
        "listing the primers in the order to report on (which "
        "can be a subset of those in the amplicon file)."
    ),
)
parser.add_argument(
    "-a",
    "--amplicons",
    nargs="+",
    metavar="TSV",
    required=True,
    help=(
        "Deduplicated bed output file(s) from isPcr. Only the "
        "first 4 columns are used."
    ),
)
parser.add_argument(
    "-o",
    "--output",
    metavar="TSV",
    required=True,
    help="Filename for tally TSV output.",
)
args = parser.parse_args()

primer_files = args.primers
fasta_files = args.fasta
bed_files = args.amplicons
tally_file = args.output


def load_primers(primer_files):
    primers = {}
    for primer_file in primer_files:
        # sys.stderr.write(f"DEBUG: Loading primer TSV file {primer_file}\n")
        for line in open(primer_file):
            if line.startswith("#") or not line.strip():
                continue
            name, fwd, rev = line.rstrip().split("\t")[:3]
            if name not in primers:
                primers[name] = {(fwd, rev)}
            elif (fwd, rev) not in primers[name]:
                primers[name].add((fwd, rev))
            else:
                sys.stderr.write(f"WARNING - duplicate line for {name}\n")
    for name in primers:
        if (count := len(primers[name])) > 1:
            sys.stderr.write(f"WARNING - cocktail of {count} pairs for {name}\n")
    return primers


primer_defs = load_primers(primer_files)
sys.stderr.write(f"Loaded {len(primer_defs)} primers\n")


def load_fasta(fasta_files):
    acc_description = {}
    for fasta_file in fasta_files:
        # sys.stderr.write(f"DEBUG: Loading reference FASTA file {fasta_file}\n")
        with open(fasta_file) as handle:
            for title, seq in SimpleFastaParser(handle):
                acc, description = title.split(None, 1)
                if acc in acc_description:
                    sys.exit(f"ERROR - Duplicate {acc} in FASTA files")
                acc_description[acc] = description
    return acc_description


if os.path.isfile(tally_file):
    sys.stderr.write(f"WARNING - Overwriting {tally_file}\n")

acc_description = load_fasta(fasta_files)
sys.stderr.write(f"Loaded lineage for {len(acc_description)} reference ids\n")

amplicon_lengths = {}
for bed_file in bed_files:
    # sys.stderr.write(f"DEBUG: Loading BED file {bed_file}\n")
    with open(bed_file) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            # We dropped column 5 (score), not checking 6 (now 5) strand etc
            acc, start, end, name, strand = line.rstrip().split("\t", 4)
            if acc not in acc_description:
                sys.exit(f"ERROR - Unexpected sequence {acc} in {bed_file}")
            # Need the lengths of the primers
            try:
                cocktail = primer_defs[name]
            except KeyError:
                sys.stderr.write(
                    f"WARNING - Ignoring unknown primer {name} in {bed_file}\n"
                )
                continue
            f_lengths = {len(f) for f, r in cocktail}
            r_lengths = {len(r) for f, r in cocktail}
            assert (
                len(f_lengths) == 1 and len(r_lengths) == 1
            ), f"Assorted lengths in {name} cocktail"
            # Do NOT add +1, the start/end are python style, len=end-start
            product_len = (
                int(end) - int(start) - list(f_lengths)[0] - list(r_lengths)[0]
            )
            # TODO - drop this and then would work even with original column 5 score?:
            assert strand in "+-", line
            try:
                amplicon_lengths[acc, name].append(product_len)
            except KeyError:
                amplicon_lengths[acc, name] = [product_len]

if not amplicon_lengths:
    sys.exit("ERROR - No amplicons loaded. Are the amplicon and primer files matched?")

with open(tally_file, "w") as handle:
    handle.write("#Sequence\tDescription\t" + "\t".join(primer_defs) + "\n")
    for acc, lineage in acc_description.items():
        fields = [acc, lineage] + [
            ";".join(str(_) for _ in sorted(amplicon_lengths.get((acc, name), [])))
            for name in primer_defs
        ]
        handle.write("\t".join(fields) + "\n")
del amplicon_lengths
sys.stderr.write(
    f"Wrote {tally_file} with {len(acc_description)} sequences vs {len(primer_defs)} primer pairs\n"
)
