#!/usr/bin/env python
#
# Based on report_len_###.py developed for an in-silico PCR
# primer evaluation project Sep 2023 to Jan 2024.
import argparse
import os
import sys

from Bio.SeqIO.FastaIO import SimpleFastaParser

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.9.1")
    sys.exit(0)

usage = """\
Parses a set of FASTA files and isPcr BED output files, to
produce a single TSV with one column per primer pair, and one
row for each taxonomic lineage inferred from the FASTA headers.
This is then used with sister script ``isPcr_lineage_report.py``
to produce tables.

Inputs:

* Primer definitions in 3-column TSV format used by isPcr
* Set of reference FASTA files where the description line is
  a semi-colon separated taxonomic lineage, as produced by
  script ``species_dedup_gbk.py`` from the source lines in NCBI
  GenBank format files.
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
    $ ./isPcr_lineage_tally.py -f refs.fasta \\
      -p primers.tsv -a amplicons.tsv -o tally.tsv

Here ``primers.tsv`` is a three-column input TSV file of primer
pair name, forward, and reverse sequences -- possibly ambiguous.
Intermediate file ``expanded.tsv`` is the expanded file of
unambiguous primers, ``amplicons.tsv`` is the isPcr output in BED
format (without column 5, score), sorted and deduplicated.
Finally ``tally.tsv`` is the output tally TSV filename.
"""

parser = argparse.ArgumentParser(
    prog="isPcr_lineage_tally.py",
    description="Produce tally of Jim Kent's isPcr results vs lineage.",
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
        "first 4 columns are used"
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
            if line.startswith("#"):
                continue
            name, fwd, rev = line.rstrip().split("\t")[:3]
            if name not in primers:
                primers[name] = (fwd, rev)
            elif primers[name] != (fwd, rev):
                sys.exit(
                    f"ERROR - inconsistent definition for {name}, "
                    f"f={primers[name][0]} r={primers[name][1]} "
                    f"versus f={fwd} r={rev}"
                )
    return primers


primer_defs = load_primers(primer_files)
sys.stderr.write(f"Loaded {len(primer_defs)} primers\n")


def load_lineages(fasta_files):
    lineage_counts = {}
    id_to_lineage = {}
    for fasta_file in fasta_files:
        # sys.stderr.write(f"DEBUG: Loading reference FASTA file {fasta_file}\n")
        with open(fasta_file) as handle:
            for title, seq in SimpleFastaParser(handle):
                del seq
                idn, lineage = title.split(None, 1)
                if idn in id_to_lineage:
                    sys.exit(f"ERROR - Duplicate {idn} in FASTA files")
                # assert lineage.startswith("Eukaryota;"), lineage
                assert ";" in lineage, lineage
                try:
                    lineage_counts[lineage] += 1
                except KeyError:
                    lineage_counts[lineage] = 1
                id_to_lineage[idn] = lineage
    return lineage_counts, id_to_lineage


if os.path.isfile(tally_file):
    sys.stderr.write(f"WARNING - Overwriting {tally_file}\n")

counts, id_to_lineage = load_lineages(fasta_files)
if len(counts) != len(id_to_lineage):
    sys.stderr.write(
        f"WARNING - Loaded {len(counts)} lineages and {len(id_to_lineage)} "
        "reference ids (expected one to one)\n"
    )
else:
    sys.stderr.write(f"Loaded lineage for {len(id_to_lineage)} reference ids\n")

amplicon_lengths = {}
for bed_file in bed_files:
    # sys.stderr.write(f"DEBUG: Loading BED file {bed_file}\n")
    with open(bed_file) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            # We dropped column 5 (score), not checking 6 (now 5) strand etc
            acc, start, end, name, strand = line.rstrip().split("\t", 4)
            try:
                fwd, rev = primer_defs[name]
            except KeyError:
                # Warning?
                continue
            # Do NOT add +1, the start/end are python style, len=end-start
            product_len = int(end) - int(start) - len(fwd) - len(rev)
            # TODO - drop this and then would work even with original column 5 score?:
            assert strand in "+-", line
            lineage = id_to_lineage[acc]
            try:
                amplicon_lengths[lineage, name].append(product_len)
            except KeyError:
                amplicon_lengths[lineage, name] = [product_len]

if not amplicon_lengths:
    sys.exit("ERROR - No amplicons loaded. Are the amplicon and primer files matched?")

with open(tally_file, "w") as handle:
    handle.write("#Lineage vs product-len\tReference\t" + "\t".join(primer_defs) + "\n")
    for lineage, count in sorted(counts.items()):
        assert count > 0, lineage
        # if all((lineage, name) not in amplicon_lengths for name in primer_defs):
        #     sys.stderr.write(f"DEBUG: No hits from {lineage}\n")
        fields = [lineage, str(count)] + [
            ";".join(str(_) for _ in sorted(amplicon_lengths.get((lineage, name), [])))
            for name in primer_defs
        ]
        handle.write("\t".join(fields) + "\n")
del amplicon_lengths
sys.stderr.write(
    f"Wrote {tally_file} with {len(counts)} lineages vs {len(primer_defs)} primer pairs\n"
)
