#!/usr/bin/env python
"""Python script for shredding contigs into fake reads.

e.g. for input into Newbler.
"""

import os
import sys
from optparse import OptionParser

from Bio import SeqIO

usage = """Basic usage: python shred_contigs.py assembly.fasta -o shredded.fasta

Multiple input FASTA files are accepted, see -h for more details.

Using Roche 454 Newbler, non-SFF input reads are limited to 1999 bp, thus
you might wish to use something like this on an Illumina assembly:

$ python shred_contigs.py other_assemby.fasta -o shredded.fasta -m 1999 -l 1999 -s 500
"""

parser = OptionParser(usage=usage)
parser.add_option(
    "-m",
    "--min-contig-len",
    dest="max_contig",
    type="int",
    help="Max contig length to reuse as is (default 2000)",
    default=2000,
)
parser.add_option(
    "-l",
    "--shred-length",
    dest="shred_length",
    type="int",
    help="Length of fake reads to generate (default 1000 bp)",
    default=1000,
)
parser.add_option(
    "-s",
    "--shred-step",
    dest="shred_step",
    type="int",
    help="Offset between fake reads (default 500 bp)",
    default=500,
)
parser.add_option(
    "-o",
    "--output",
    dest="output_filename",
    help="FASTA output filename for fake reads (required)",
    default=None,
    metavar="FILE",
)
(options, args) = parser.parse_args()
if not args:
    sys.exit("Requires at least one input FASTA filename\n\n" + usage)

max_contig = int(options.max_contig)
shred_length = int(options.shred_length)
shred_step = int(options.shred_step)
output_fasta = options.output_filename

if shred_step < 1:
    sys.exit("Shred step should be positive")
if shred_length < shred_step:
    sys.exit("Shred step should be less than shred length")

print("Accepting contigs up to length %i as they are (option -m)" % max_contig)
print(
    "Shredding longer contigs into reads of %i bp (option -l), step %i (option -s)"
    % (shred_length, shred_step)
)

for assembly_fasta in args:
    if not os.path.isfile(assembly_fasta):
        sys.exit("Assembly FASTA file not found: %r" % assembly_fasta)


def shred(input_filename):
    global as_is, shredded
    for record in SeqIO.parse(input_filename, "fasta"):
        length = len(record)
        if length <= max_contig:
            as_is += 1
            yield record
        else:
            # Shred it!
            shredded += 1
            for i, start in enumerate(range(0, length - shred_step, shred_step)):
                fragment = record[start : start + shred_length]
                fragment.id = "%s_fragment%i" % (record.id, i + 1)
                yield fragment


count = 0
as_is = 0
shredded = 0
with open(output_fasta, "w") as output_handle:
    for assembly_fasta in args:
        count += SeqIO.write(shred(assembly_fasta), output_handle, "fasta")
print("Shredded %i FASTA files, containing %i contigs" % (len(args), as_is + shredded))
print(
    "Accepted %i short contigs, shredded %i long contigs, giving %i reads"
    % (as_is, shredded, count)
)
