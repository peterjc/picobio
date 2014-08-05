#!/usr/bin/env python
"""Python script for shredding contigs into fake reads.

e.g. for input into Newbler.
"""

import os
import sys
import warnings
import tempfile
import shutil
from optparse import OptionParser
from Bio import SeqIO

usage = """Basic usage: shred_contigs.py assembly.fasta -o dedup_output.fasta
"""

def stop_err(msg, error_level=1):
    """Print error message to stdout and quit with given error level."""
    sys.stderr.write("%s\n" % msg.rstrip())
    sys.exit(error_level)

parser = OptionParser(usage=usage)
parser.add_option("-m", "--min-contig-len", dest="max_contig", type="int",
                  help="Max contig length to re-use as is (default 2000)",
                  default=2000)
parser.add_option("-l", "--shred-length", dest="shred_length", type="int",
                  help="Length of fake reads to generate (default 1000 bp)",
                  default=1000)
parser.add_option("-s", "--shred-step", dest="shred_step", type="int",
                  help="Offset between fake reads (default 500 bp)",
                  default=500)
parser.add_option("-o", "--output", dest="output_filename",
                  help="FASTA output filename for fake reads (required)",
                  default=None,
                  metavar="FILE")
(options, args) = parser.parse_args()
if not args:
    stop_err("Requires at least one input FASTA filename\n\n" + usage)

max_contig = int(options.max_contig)
shred_length = int(options.shred_length)
shred_step = int(options.shred_step)
output_fasta = options.output_filename

if shred_step < 1:
    stop_err("Shred step should be positive")

for assembly_fasta in args:
    if not os.path.isfile(assembly_fasta):
        stop_err("Assembly FASTA file not found: %r" % assembly_fasta)

def shred(input_filename):
    for record in SeqIO.parse(input_filename, "fasta"):
        length = len(record)
        if length <= max_contig:
            yield record
        else:
            #Shred it!
            for i, start in enumerate(range(0, length, shred_step)):
                fragment = record[start:start+shred_length]
                fragment.id = "%s_fragment%i" % (record.id, i+1)
                yield fragment


count = 0
with open(output_fasta, "w") as output_handle:
    for assembly_fasta in args:
        count += SeqIO.write(shred(assembly_fasta), output_handle, "fasta")
print("Shreded %i FASTA files, giving %i reads" % (len(args), count))
