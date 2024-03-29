#!/usr/bin/env python
#
# Based on tally_coliphage_primers_###.py developed
# while evaluating and designing primers for coliphage
# in August 2023 for a RESAS project.
import argparse
import sys

import numpy as np
import pandas as pd
import seaborn as sns
from Bio.SeqIO.FastaIO import SimpleFastaParser
from matplotlib import pyplot as plt

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)

usage = """\
This script takes the output of Jim Kent's isPcr tool in FASTA format as
input. It computes a matrix with accessions as rows, and primer sets as
columns, with amplicon length as the value (or zero). This is output as a TSV
file, and plotted as a heatmap output as a PDF file.

Conflicting amplicon lengths (which requires multiple products from a given
reference sequence) is considered to be an error, and such primers are
discarded.
"""

parser = argparse.ArgumentParser(
    prog="plot_ispcr.py",
    description="Plot heatmap from output of Jim Kent's isPcr tool.",
    epilog=usage,
)
parser.add_argument(
    "-i",
    "--input",
    dest="input",
    nargs="+",
    default=["/dev/stdin"],
    metavar="FILE",
    help="Input FASTA filename(s) from isPcr, default stdin.",
)
parser.add_argument(
    "-t",
    "--target",
    type=int,
    default="0",
    metavar="INTEGER",
    help="Target amplicon length (for a default divergent color scheme).",
)
parser.add_argument(
    "-m",
    "--mincount",
    type=int,
    default="5",
    metavar="INTEGER",
    help="Minimum number of hits required to report on a primer pair (default 5).",
)
parser.add_argument(
    "-x",
    "--maxlength",
    type=int,
    default="1000",
    metavar="INTEGER",
    help="Maximum amplicon size to report on a primer pair (default 1000).",
)
parser.add_argument(
    "--labelallprimers",
    action="store_true",
    default=False,
    help="Force captions for all primers (may overlap).",
)
parser.add_argument(
    "--labelallseqs",
    action="store_true",
    default=False,
    help="Force captions for all reference sequences (may overlap).",
)
parser.add_argument(
    "--vmin",
    type=int,
    default="-1",
    metavar="INTEGER",
    help="Minimum amplicon length for colour scheme (default automatic)",
)
parser.add_argument(
    "--vmax",
    type=int,
    default="-1",
    metavar="INTEGER",
    help="Maximum amplicon length for colour scheme (default automatic)",
)
# parser.add_argument(
#    "-d",
#    "--database",
#    nargs="+",
#    metavar="FILE",
#    help="FASTA filename(s) given to isPcr, or a text file containing those filenames.",
# )
parser.add_argument(
    "-o",
    "--output",
    dest="output",
    required=True,
    metavar="STEM",
    help="Output filename stem (required).",
)

if len(sys.argv) == 1:
    sys.exit("ERROR: Invalid command line, try -h or --help.")
options = parser.parse_args()
if options.target and options.target < 40:
    sys.exit(f"ERROR: If used target length should be over 40, not {options.target}")
if options.maxlength and options.maxlength < 40:
    sys.exit(f"ERROR: If used max-length should be over 40, not {options.maxlength}")
if options.target > options.maxlength:
    sys.exit(f"ERROR: Max length {options.maxlength} < target length {options.target}")


def main(
    pcr_results,
    output_stem,
    min_count=0,
    max_length=1000,
    target_length=0,
    label_all_primers=False,
    label_all_seqs=False,
    vmin=None,
    vmax=None,
):
    products = {}
    rejected = set()
    for pcr_file in pcr_results:
        sys.stderr.write(f"Loading isPcr predictions from {pcr_file}\n")
        with open(pcr_file) as handle:
            for title, seq in SimpleFastaParser(handle):
                acc_loc, primer_name, amplicon_size, forward, reverse, acc_name = (
                    title + " "
                ).split(" ", 5)
                assert amplicon_size.endswith("bp"), title
                amplicon_size = int(amplicon_size[:-2])
                acc = acc_loc.split(":", 1)[0]
                if primer_name not in rejected and amplicon_size > max_length:
                    sys.stderr.write(
                        f"Rejecting {primer_name} as can give too long an amplicon.\n"
                    )
                    rejected.add(primer_name)
                if primer_name not in rejected and amplicon_size != products.get(
                    (acc, primer_name), amplicon_size
                ):
                    sys.stderr.write(
                        f"Rejecting {primer_name} for multiple product sizes.\n"
                    )
                    rejected.add(primer_name)
                products[acc, primer_name] = amplicon_size
    if not products:
        sys.exit("ERROR: Loaded no in-silico PCR products")
    sys.stderr.write(
        f"Loaded {len(products)} in-silico PCR results, "
        f"max amplicon size {max(products.values())}\n"
    )

    hits = sorted({acc for acc, primer_name in products})
    sys.stderr.write(f"Was able to amplify {len(hits)} accessions in all\n")

    primers = sorted({primer_name for acc, primer_name in products})
    sys.stderr.write(f"Have amplification from {len(primers)} primers\n")

    if rejected:
        sys.stderr.write(
            f"Rejecting {len(rejected)} primers due to multiple or large product sizes\n"
        )
    for primer_name in primers:
        count = sum(1 for acc in hits if (acc, primer_name) in products)
        if primer_name in rejected or count < min_count:
            if primer_name not in rejected:
                sys.stderr.write(f"Culling {primer_name} as {count} is too few hits\n")
            for acc in hits:
                if (acc, primer_name) in products:
                    del products[acc, primer_name]
    if not products:
        sys.exit("ERROR: No primer/accession pairs accepted")

    # Updates hits list as after dropping primers some accessions may have no hits:
    hits = sorted({acc for acc, primer_name in products})
    primers = sorted({primer_name for acc, primer_name in products})
    sys.stderr.write(
        f"Now have {len(primers)} primers vs {len(hits)} accessions, "
        f"max amplicon size {max(products.values())}\n"
    )

    # Tabular output - todo, apply the same ordering as the heatmap?
    with open(output_stem + ".tsv", "w") as out_handle:
        out_handle.write("#\t" + "\t".join(primers) + "\n")
        # for acc, row in zip(hits, as_array):
        #    out_handle.write(acc + "\t" + "\t".join(str(_) for _ in row) + "\n")
        for acc in hits:
            out_handle.write(
                acc
                + "\t"
                + "\t".join(
                    str(products.get((acc, primer_name), 0)) for primer_name in primers
                )
                + "\n"
            )
    sys.stderr.write(f"Wrote {output_stem}.tsv\n")

    as_array = np.array(
        [
            [products.get((acc, primer_name), 0) for primer_name in primers]
            for acc in hits
        ],
        np.int16,
    )
    data_frame = pd.DataFrame(as_array, columns=primers, index=hits)

    cluster_grid = sns.clustermap(
        data_frame,  # as_array,
        mask=(as_array == 0),
        # I like the default colour map in target mode, otherwise "flare":
        center=target_length if target_length else None,
        cmap=None if target_length else "flare",
        row_cluster=(len(hits) > 1),
        col_cluster=(len(primers) > 1),
        xticklabels=True if label_all_primers else "auto",  # may overlap
        yticklabels=True if label_all_seqs else "auto",  # may overlap
        vmin=vmin,
        vmax=vmax,
        robust=True,  # ignored if both vmin and vmax is used
    )
    # Does this work for smaller font?:
    plt.setp(cluster_grid.ax_heatmap.get_xticklabels(), fontsize=8)  # For x axis

    cluster_grid.savefig(output_stem + ".pdf")
    sys.stderr.write(f"Wrote {output_stem}.pdf\n")
    cluster_grid.savefig(output_stem + ".png")
    sys.stderr.write(f"Wrote {output_stem}.png\n")


main(
    options.input,
    options.output,
    min_count=options.mincount,
    max_length=options.maxlength,
    target_length=options.target,
    label_all_primers=options.labelallprimers,
    label_all_seqs=options.labelallseqs,
    vmin=options.vmin if options.vmin != -1 else None,
    vmax=options.vmax if options.vmax != -1 else None,
)
